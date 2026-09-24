"""Stand-in for pythonnet: records how .NET was asked to load.
(The knobs: ThermoFisher/CommonCore/RawFileReader/__init__.py.)"""
import json
import os

_LOADED = []


def load(runtime=None, **params):
    log = os.environ.get("FAKE_PYTHONNET_LOG")
    if log:
        with open(log, "a") as fh:
            fh.write(json.dumps({"runtime": runtime, "params": params,
                                 "DOTNET_ROOT": os.environ.get("DOTNET_ROOT")}) + "\n")
    if os.environ.get("FAKE_PYTHONNET_LOAD_FAIL"):
        raise RuntimeError(f"Failed to create a .NET runtime ({runtime}) using the parameters "
                           f"{params}.")
    root = params.get("dotnet_root")
    if root and not os.path.isdir(os.path.join(root, "host", "fxr")):
        raise RuntimeError(f"Could not find a suitable hostfxr library in {root}")
    _LOADED.append(runtime)


def get_runtime_info():
    return None
