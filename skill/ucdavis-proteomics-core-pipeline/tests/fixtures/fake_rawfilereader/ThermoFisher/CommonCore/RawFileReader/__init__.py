"""
Stand-ins for pythonnet, `clr` and the `ThermoFisher.CommonCore.*` .NET namespaces, so
`scripts/thermo_resolution.py` runs for real in tests without .NET or Thermo's DLLs. Put this
folder on PYTHONPATH (and run Python with `-S`, so an installed pythonnet cannot win).

A fake .raw is described by `<raw>.scans.json` beside it:
  {"n": 5000, "cycle": [{"filter": "FTMS + p NSI Full ms [350-1500]",
                         "trailer": {"Orbitrap Resolution:": "60000"}}, ...]}
scan i (1-based) is cycle[(i - 1) % len(cycle)]; {"error": "<message>"} makes it unopenable.

Environment knobs:
  FAKE_PYTHONNET_LOG=<file>     one JSON line per pythonnet.load(): runtime, params, DOTNET_ROOT
  FAKE_PYTHONNET_LOAD_FAIL=1    load() raises as pythonnet does when no runtime can be created
  FAKE_RFR_LOG=<file>           one line per RawFileReaderAdapter.FileFactory(path)
  FAKE_RFR_SLEEP_ON=<names>     FileFactory sleeps FAKE_RFR_SLEEP seconds (default 30) on a
                                path that has one of these (comma-separated) as a component
  FAKE_RFR_CRASH_ON=<names>     ... kills the process (exit 134) on such a path: a native crash
"""
import json
import os
import time


class _Text:
    def __init__(self, text):
        self._text = text

    def ToString(self):
        return self._text


class _Trailer:
    def __init__(self, trailer):
        self.Labels = list(trailer)
        self.Values = [str(v) for v in trailer.values()]
        self.Length = len(self.Labels)


class _Header:
    def __init__(self, n):
        self.FirstSpectrum, self.LastSpectrum = 1, n


class _Error:
    def __init__(self, message):
        self.ErrorMessage = message


class _Raw:
    def __init__(self, spec):
        self.IsError = "error" in spec
        self.FileError = _Error(spec.get("error", ""))
        self._cycle = spec.get("cycle") or []
        self.RunHeaderEx = _Header(int(spec.get("n", len(self._cycle))))

    def SelectInstrument(self, device, index):
        pass

    def _scan(self, scan):
        return self._cycle[(scan - 1) % len(self._cycle)]

    def GetFilterForScanNumber(self, scan):
        return _Text(self._scan(scan)["filter"])

    def GetTrailerExtraInformation(self, scan):
        return _Trailer(self._scan(scan).get("trailer", {}))

    def Dispose(self):
        pass


def _hit(knob, path):
    names = [n for n in os.environ.get(knob, "").split(",") if n]
    return any(n in path.split(os.sep) for n in names)


class RawFileReaderAdapter:
    @staticmethod
    def FileFactory(path):
        log = os.environ.get("FAKE_RFR_LOG")
        if log:
            with open(log, "a") as fh:
                fh.write(path + "\n")
        if _hit("FAKE_RFR_CRASH_ON", path):
            os._exit(134)
        if _hit("FAKE_RFR_SLEEP_ON", path):
            time.sleep(float(os.environ.get("FAKE_RFR_SLEEP", "30")))
        try:
            with open(path + ".scans.json") as fh:
                spec = json.load(fh)
        except (OSError, ValueError) as e:
            spec = {"error": f"no fake scan spec: {e}"}
        return _Raw(spec)
