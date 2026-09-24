"""Stand-in for pythonnet's `clr`: AddReference finds <name>.dll on sys.path.
(The knobs: ThermoFisher/CommonCore/RawFileReader/__init__.py.)"""
import os
import sys


class _Version:
    def ToString(self):
        return "8.0.6.0"


class _Name:
    Version = _Version()


class _Assembly:
    def GetName(self):
        return _Name()


def AddReference(name):
    for d in sys.path:
        if d and os.path.isfile(os.path.join(d, name + ".dll")):
            return _Assembly()
    raise FileNotFoundError(f"Unable to find assembly '{name}'.")
