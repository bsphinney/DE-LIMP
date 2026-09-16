#!/usr/bin/env python3
"""
A stand-in for ThermoRawFileParser that replays output captured from the real one.

It is deliberately strict about the command line, because the defect it exists to catch was
a command line: detect_acquisition.py called `ThermoRawFileParser metadata -i <raw>` and
`ThermoRawFileParser query -i <raw>`, and the real parser (v.2.0.0-dev, 2.0.0.0, on HIVE,
srun job 23510571) rejected both, with exit 255 and these first lines on stderr:

    Unexpected extra arguments                     <- there is no `metadata` subcommand
    -s, --scans: specify a valid scan range        <- `query` requires -n

The grammar reproduced here (short options, -x=value) is accepted by every release from
v1.4.0 through v.2.0.0-dev (MainClass.cs; only the long option names moved, in 1.4.4):
  metadata only:  -i=<raw> -m=0 -o=<dir>    -> <dir>/<stem>-metadata.json
  spectra query:  query -i=<raw> -n=<scans> -b=<file>  -> <dir of file>/<stem of file>.json
  version:        --version                  -> "2.0.0.0"

Replay: the raw's file stem picks the fixture (<stem>-metadata.json / <stem>.query.json in
this directory). The raw must exist, as the real parser checks File.Exists first.

Test knobs (environment):
  FAKE_TRFP_LOG=<file>      append each argv as one JSON line, so tests can assert the call
  FAKE_TRFP_FAIL=<mode>     make `metadata` or `query` fail the way a damaged raw does:
                            ERROR line on STDOUT (log4net's console appender), exit 1
  FAKE_TRFP_GARBAGE=query   write a query file that is not JSON (truncated output)
  FAKE_TRFP_NO_FILTER=1     drop the filter-string attribute (MS:1000512), as 1.4.0-1.4.4
                            write it: those builds report isolation windows but not the
                            data-dependent flag
"""
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
VALUE_OPTS = {"i": "input", "n": "scans", "b": "output", "o": "output_directory",
              "m": "metadata", "f": "format", "l": "logging", "c": "metadata_output_file"}
LONG_TO_SHORT = {v: k for k, v in VALUE_OPTS.items()}


def usage_error(first_line):
    sys.stderr.write(first_line + "\nError - usage is:\n  -h, --help   Prints out the options.\n")
    sys.exit(255)


def parse(args):
    opts, extra = {}, []
    i = 0
    while i < len(args):
        a = args[i]
        m = re.match(r"^--?([A-Za-z_]+)(?:[=:](.*))?$", a, re.S)
        if not m:
            extra.append(a)
            i += 1
            continue
        name, val = m.group(1), m.group(2)
        name = LONG_TO_SHORT.get(name, name)
        if name in VALUE_OPTS and val is None:
            i += 1
            val = args[i] if i < len(args) else None
        opts[name] = True if val is None else val
        i += 1
    return opts, extra


def main():
    argv = sys.argv[1:]
    log = os.environ.get("FAKE_TRFP_LOG")
    if log:
        with open(log, "a") as fh:
            fh.write(json.dumps(argv) + "\n")

    if argv in (["--version"], ["-v"]):
        print("2.0.0.0")
        return 0

    sub = argv[0] if argv and argv[0] in ("query", "xic") else None
    opts, extra = parse(argv[1:] if sub else argv)
    if extra:
        usage_error("Unexpected extra arguments" if not sub else "unexpected extra arguments")

    raw = opts.get("i")
    if not raw or raw is True:
        usage_error("-i, --input: specify an input file")
    if not os.path.isfile(raw):
        usage_error("-i, --input: specify a valid RAW file location")
    stem = os.path.splitext(os.path.basename(raw))[0]
    fail = os.environ.get("FAKE_TRFP_FAIL", "")

    if sub == "query":
        scans = opts.get("n")
        if not scans or scans is True or not re.match(r"^[\d,\-\s]+$", scans):
            usage_error("-s, --scans: specify a valid scan range")
        if fail == "query":
            print("2026-09-16 16:09:00 ERROR RAW file cannot be processed because of an error - "
                  "file is corrupt")
            return 1
        out = opts.get("b")
        if out and out is not True:
            dest = os.path.join(os.path.dirname(os.path.abspath(out)),
                                os.path.splitext(os.path.basename(out))[0] + ".json")
        else:
            dest = os.path.join(os.path.dirname(os.path.abspath(raw)), stem + ".json")
        with open(dest, "w") as fh:
            if os.environ.get("FAKE_TRFP_GARBAGE") == "query":
                fh.write('[{"mzs":[350.0,')
            else:
                with open(os.path.join(HERE, stem + ".query.json")) as src:
                    spectra = json.load(src)
                if os.environ.get("FAKE_TRFP_NO_FILTER"):
                    for spec in spectra:
                        spec["attributes"] = [a for a in spec["attributes"]
                                              if a.get("accession") != "MS:1000512"]
                json.dump(spectra, fh)
        return 0

    if sub is None and opts.get("m") is not None:
        if opts["m"] not in ("0", "json", "JSON"):
            usage_error("-m, --metadata: this stand-in only replays JSON metadata")
        if fail == "metadata":
            print("2026-09-16 16:09:00 ERROR RAW file cannot be processed because of an error - "
                  "file is corrupt")
            return 1
        outdir = opts.get("o")
        if not outdir or outdir is True:
            outdir = os.path.dirname(os.path.abspath(raw))
        if not os.path.isdir(outdir):
            usage_error("-o, --output: specify a valid output directory")
        with open(os.path.join(HERE, stem + "-metadata.json")) as src, \
                open(os.path.join(outdir, stem + "-metadata.json"), "w") as fh:
            fh.write(src.read())
        return 0

    usage_error("this stand-in replays only metadata (-m=0) and query (-n) calls")


if __name__ == "__main__":
    sys.exit(main())
