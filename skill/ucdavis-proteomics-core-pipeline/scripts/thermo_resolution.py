#!/usr/bin/env python3
"""
thermo_resolution.py  --  the Orbitrap MS1/MS2 RESOLUTION of Thermo .raw files, read from the
scan trailer with Thermo's own RawFileReader (.NET), through pythonnet.

Why this exists: estimate_params.py pins DIA-NN's documented Orbitrap mass tolerances from the
MS1/MS2 resolution (240k->4, 120k->7, 60k->10, 30k->15 ppm). Without it an Orbitrap classes as
`orbitrap_generic`: DIA-NN calibrates per run and the 5-step parallel chain declines the cfg
(gabrig 2026-09-23, 15 Fusion Lumos .raw). The instrument writes the resolution into every
scan's trailer, but ThermoRawFileParser never outputs it -- TRFP master cf548e4 reads trailer
keys internally and writes no resolution to metadata, mzML, MGF or `query`; its
"mass resolution" (MS:1000011) is a generic 0.5. Verified on HIVE: a Fusion Lumos writes
'Orbitrap Resolution:' (jobs 23989081, 23989291: MS1 60000, MS2 15000); an Exploris 480 writes
'FT Resolution:' at both levels (job 23989339: two HeLa DIA runs 120000/15000, a DDA sample
60000/15000). The instrument method TEXT is not used, only because it is not always there: the
Exploris 480 blank first probed (Exploris480/Blanx/Ex1900V600nL4B-...raw, 1,973 scans) has
InstrumentMethodsCount 0. Where it is there it agrees -- the three Exploris runs of job 23989339
have 2 methods each, whose "Orbitrap Resolution = " lines match the trailer exactly (and use a
different label from the trailer's 'FT Resolution:', so a port must accept both).

detect_acquisition.py runs this as a SUBPROCESS -- loading CoreCLR into a process, or a native
crash inside it, must not be able to take acquisition detection down -- and passes the
DOTNET_ROOT it resolved (a root with Microsoft.NETCore.App 8: the RawFileReader DLLs target
.NETCoreApp 8.0). The DLLs are the ones ThermoRawFileParser ships beside its executable.

Usage:
  thermo_resolution.py --dll-dir DIR RAW [RAW ...]   one JSON object per RAW, one per line:
      {"file", "ms1_resolution", "ms2_resolution", "ms2_analyzer", "ms1_key", "ms2_key",
       "ms1_scan", "ms2_scan", "ms1_filter", "ms2_filter", "scans_walked", "scan_types",
       "reader", "note"}
      A resolution is an int, or null + the why in note. ms2_analyzer is "FTMS", "ITMS" (an
      ion-trap MS2: no Orbitrap resolution exists for it), "mixed" (both in one run -- a
      Tribrid decision-tree method, HCD-OT for high charge states and CID-IT for 2+; the MS2
      resolution then covers the Orbitrap MS2 only) or null. ms1_filter/ms2_filter are
      the full filter strings of the scans the values came from, and scan_types counts every
      scan type in the walk (filter up to its mass range, precursor m/z dropped): the evidence
      that the values are the method's own MS1 and MS2, not a SIM, a second survey or another
      scan type.
  thermo_resolution.py --dll-dir DIR --check         {"ok", "reader", "note"}: can it load at all
Never raises. Exit 0 when every RAW was answered (answers may be null), 1 when the runtime or the
DLLs could not be loaded (every RAW still gets its line, with the reason), 2 on bad usage.
"""
import json
import os
import re
import sys
import tempfile
import time

DLLS = ("ThermoFisher.CommonCore.RawFileReader", "ThermoFisher.CommonCore.Data")
# Trailer labels as the instruments write them, colon included.
RESOLUTION_KEYS = ("Orbitrap Resolution:", "FT Resolution:")
# The resolution is a method setting, the same on every scan of a level, so a few cycles is
# enough -- but a DDA run's first minutes (void volume) can be MS1 only, so the walk starts
# mid-run, where the run is at steady state (as detect_acquisition's TRFP query does).
WALK_SCANS = 400
WALK_SECONDS = 20             # per file (~3 s measured); inside detect_acquisition's 30 s/file
PYTHONNET_MISSING = (
    "pythonnet is not installed for this Python, so the Orbitrap resolution could not be read "
    "from the scan trailer. `bash scripts/setup.sh` installs it into the pipeline env "
    "(conda-forge `pythonnet`), or `pip install pythonnet`; then run step 2 with that env's "
    "python3 (`source ~/.proteomics-pipeline/activate.sh`)")


def _emit(obj):
    print(json.dumps(obj), flush=True)       # flushed: a later crash must not lose this line


def scan_level(filter_text):
    """(analyzer, MS order, is a Full scan) from a Thermo scan filter string.

    "FTMS + p NSI Full ms [350.0000-1500.0000]"             -> ("FTMS", 1, True)
    "FTMS + p NSI Full ms2 368.0000@hcd35.00 [250.0-746.0]" -> ("FTMS", 2, True)
    "ITMS + c NSI r d Full ms2 572.32@cid35.00 [...]"       -> ("ITMS", 2, True)
    "FTMS + p NSI Full msx ms2 ..." (multiplexed)           -> ("FTMS", 2, True)
    """
    head = (filter_text or "").split("[", 1)[0].split()
    analyzer = head[0].upper() if head else ""
    order = None
    for tok in head:
        m = re.fullmatch(r"ms(\d*)", tok)
        if m:
            order = int(m.group(1) or 1)
    return analyzer, order, "Full" in head


def scan_type(filter_text):
    """The filter string's scan TYPE: up to its mass range, each precursor m/z dropped.

    "FTMS + p NSI Full ms2 368.0000@hcd35.00 [250.0000-746.0000]"
        -> "FTMS + p NSI Full ms2 @hcd35.00"
    """
    head = (filter_text or "").split("[", 1)[0].split()
    return " ".join(("@" + t.split("@", 1)[1]) if "@" in t else t for t in head)


def trailer_resolution(labels, values):
    """(key, resolution) from one scan's trailer, or (None, None)."""
    for label, value in zip(labels, values):
        label = str(label).strip()
        if label in RESOLUTION_KEYS:
            try:
                res = int(float(str(value).strip()))
            except ValueError:
                continue
            if res > 0:
                return label, res
    return None, None


def _runtime_config(root):
    """(runtimeconfig.json path, version) for the highest Microsoft.NETCore.App in `root`, or
    (None, None).

    Written here rather than left to clr_loader: its own discovery takes the runtime list from
    `dotnet --list-runtimes` of whatever dotnet is on PATH even when DOTNET_ROOT names another
    root, and would pair one root's hostfxr with another root's runtime.
    """
    shared = os.path.join(root, "shared", "Microsoft.NETCore.App")
    try:
        versions = [v for v in os.listdir(shared) if re.fullmatch(r"\d+\.\d+\.\d+", v)]
    except OSError:
        return None, None
    if not versions:
        return None, None
    best = max(versions, key=lambda v: tuple(int(x) for x in v.split(".")))
    fd, path = tempfile.mkstemp(prefix="thermo_resolution_", suffix=".runtimeconfig.json")
    with os.fdopen(fd, "w") as fh:
        json.dump({"runtimeOptions": {"tfm": "net" + ".".join(best.split(".")[:2]),
                                      "framework": {"name": "Microsoft.NETCore.App",
                                                    "version": best}}}, fh)
    return path, best


def load_reader(dll_dir):
    """(RawFileReaderAdapter, Device, reader description, None) or (None, None, None, why)."""
    try:
        import pythonnet
    except ImportError:
        return None, None, None, PYTHONNET_MISSING
    missing = [d for d in DLLS if not os.path.isfile(os.path.join(dll_dir, d + ".dll"))]
    if missing:
        return None, None, None, (f"no {', '.join(d + '.dll' for d in missing)} in {dll_dir}")
    root = os.environ.get("DOTNET_ROOT")
    params, net, cfg = {}, "?", None
    if root:
        cfg, net = _runtime_config(root)
        if not cfg:
            return None, None, None, (f"DOTNET_ROOT={root} has no Microsoft.NETCore.App runtime "
                                      "(shared/Microsoft.NETCore.App/<version>)")
        params = {"runtime_config": cfg, "dotnet_root": root}
    try:
        try:
            pythonnet.load("coreclr", **params)
        finally:
            # hostfxr has read it by the time load() returns (or fails): no temp file left behind
            if cfg:
                try:
                    os.remove(cfg)
                except OSError:
                    pass
        import clr
        sys.path.append(dll_dir)
        versions = []
        for name in DLLS:
            asm = clr.AddReference(name)
            versions.append(str(asm.GetName().Version.ToString()))
        from ThermoFisher.CommonCore.RawFileReader import RawFileReaderAdapter
        from ThermoFisher.CommonCore.Data.Business import Device
    except Exception as e:              # the runtime, the DLLs: anything .NET can throw
        return None, None, None, (f"pythonnet could not load .NET (CoreCLR"
                                  f"{', DOTNET_ROOT=' + root if root else ''}) or the "
                                  f"RawFileReader DLLs from {dll_dir}: "
                                  f"{type(e).__name__}: {str(e).strip()[:300]}")
    try:
        from importlib.metadata import version as _pkg_version
        pn = _pkg_version("pythonnet")
    except Exception:
        pn = "?"
    reader = (f"ThermoFisher.CommonCore.RawFileReader {versions[0]} via pythonnet {pn} "
              f"(.NET {net}"
              f"{', DOTNET_ROOT=' + root if root else ''}) from {dll_dir}")
    return RawFileReaderAdapter, Device, reader, None


def blank(path, note=None):
    """The per-file object with nothing read yet (the keys: see the module docstring)."""
    return {"file": path, "ms1_resolution": None, "ms2_resolution": None, "ms2_analyzer": None,
            "ms1_key": None, "ms2_key": None, "ms1_scan": None, "ms2_scan": None,
            "ms1_filter": None, "ms2_filter": None, "scans_walked": None, "scan_types": {},
            "reader": None, "note": note}


def read_one(adapter, device, path):
    """The per-file object (see the module docstring). Never raises."""
    out = blank(path)
    raw = None
    try:
        raw = adapter.FileFactory(path)
        if getattr(raw, "IsError", False):
            out["note"] = f"RawFileReader could not open it: {raw.FileError.ErrorMessage}"
            return out
        raw.SelectInstrument(device.MS, 1)
        first, last = int(raw.RunHeaderEx.FirstSpectrum), int(raw.RunHeaderEx.LastSpectrum)
        start = max(first, (first + last) // 2 - WALK_SCANS // 2)
        stop = min(last, start + WALK_SCANS - 1)
        seen = {1: {}, 2: {}}            # level -> {resolution: (key, first scan, filter)}
        types = out["scan_types"]
        itms_ms2 = None
        n_ms2 = {"FTMS": 0, "ITMS": 0}     # MS2 scans by analyzer, whatever their trailer says
        t0 = time.monotonic()
        scan = start
        for scan in range(start, stop + 1):
            if time.monotonic() - t0 > WALK_SECONDS:
                out["note"] = (f"stopped after {WALK_SECONDS} s at scan {scan} of the "
                               f"{start}-{stop} walk (slow storage?)")
                break
            filt = str(raw.GetFilterForScanNumber(scan).ToString())
            types[scan_type(filt)] = types.get(scan_type(filt), 0) + 1
            analyzer, order, full = scan_level(filt)
            if order == 2 and analyzer in n_ms2:
                n_ms2[analyzer] += 1
            if order == 2 and analyzer == "ITMS":
                itms_ms2 = itms_ms2 or scan
                continue
            if analyzer != "FTMS" or order not in (1, 2) or (order == 1 and not full):
                continue
            tr = raw.GetTrailerExtraInformation(scan)
            key, res = trailer_resolution(list(tr.Labels), list(tr.Values))
            if res is not None:
                seen[order].setdefault(res, (key, scan, filt))
        out["scans_walked"] = f"{start}-{scan}"
        # BOTH analyzers is not "FTMS": searched as Orbitrap MS2, the ion-trap spectra would be
        # matched at +/-10 ppm and lost without a word.
        out["ms2_analyzer"] = ("mixed" if n_ms2["FTMS"] and n_ms2["ITMS"] else
                               "FTMS" if n_ms2["FTMS"] or seen[2] else
                               "ITMS" if n_ms2["ITMS"] else None)
        notes = [out["note"]] if out["note"] else []
        if out["ms2_analyzer"] == "mixed":
            ms2_types = ", ".join(f"{n} x '{t}'" for t, n in sorted(types.items())
                                  if scan_level(t)[1] == 2)
            notes.append(f"MS2 is MIXED: {n_ms2['FTMS']} Orbitrap (FTMS) and {n_ms2['ITMS']} "
                         f"ion-trap (ITMS) MS2 scans in scans {start}-{scan} ({ms2_types}); "
                         f"the MS2 resolution covers the Orbitrap MS2 only")
        for order, tag in ((1, "ms1"), (2, "ms2")):
            if len(seen[order]) == 1:
                res, (key, first_scan, filt) = next(iter(seen[order].items()))
                out.update({f"{tag}_resolution": res, f"{tag}_key": key,
                            f"{tag}_scan": first_scan, f"{tag}_filter": filt})
            elif len(seen[order]) > 1:
                # A method with two resolutions at one level: no single value is "the" one.
                notes.append(f"MS{order} scans {start}-{scan} carry several resolutions ("
                             + ", ".join(f"{r} from scan {s}" for r, (_k, s, _f)
                                         in sorted(seen[order].items())) + ")")
            elif order == 2 and itms_ms2 and not n_ms2["FTMS"]:
                notes.append(f"the MS2 scans are ion trap (ITMS, first at scan {itms_ms2}): "
                             f"no Orbitrap resolution, and DIA-NN's Orbitrap table does not "
                             f"apply to them")
            else:
                notes.append(f"no FTMS {'Full ms' if order == 1 else 'ms2'} scan with a "
                             f"resolution in the trailer ({' / '.join(RESOLUTION_KEYS)}) in "
                             f"scans {start}-{scan}")
        out["note"] = "; ".join(notes) or None
    except Exception as e:
        out["note"] = f"reading the scan trailers failed: {type(e).__name__}: {str(e)[:300]}"
    finally:
        try:
            if raw is not None:
                raw.Dispose()
        except Exception:
            pass
    return out


def main(argv):
    usage = ("usage: thermo_resolution.py --dll-dir DIR RAW [RAW ...]\n"
             "       thermo_resolution.py --dll-dir DIR --check")
    if len(argv) < 3 or argv[0] != "--dll-dir":
        print(usage, file=sys.stderr)
        return 2
    dll_dir, rest = argv[1], argv[2:]
    check = rest == ["--check"]
    adapter, device, reader, why = load_reader(dll_dir)
    if check:
        _emit({"ok": why is None, "reader": reader, "note": why})
        return 0 if why is None else 1
    for path in rest:
        if why:
            _emit(blank(path, why))
            continue
        rec = read_one(adapter, device, path)
        rec["reader"] = reader
        _emit(rec)
    return 0 if why is None else 1


if __name__ == "__main__":
    try:
        sys.exit(main(sys.argv[1:]))
    except Exception as e:                # never a traceback instead of an answer
        print(f"thermo_resolution.py: {type(e).__name__}: {e}", file=sys.stderr)
        sys.exit(1)
