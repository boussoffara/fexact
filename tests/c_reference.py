"""Loads the original C extension (see legacy/setup.py) for comparisons.

The C module is also called `fexact`, so it is loaded from its file and kept
out of `sys.modules`, next to the Rust one. `FEXACT_C` may name the file.
"""
import glob
import importlib.util
import os
import sys

LEGACY = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "legacy")


def load():
    """The C `fexact` module, or None if it hasn't been built."""
    path = os.environ.get("FEXACT_C")
    if not path:
        built = glob.glob(os.path.join(LEGACY, "fexact.*.so")) + glob.glob(
            os.path.join(LEGACY, "fexact.*.pyd")
        )
        if not built:
            return None
        path = built[0]
    spec = importlib.util.spec_from_file_location("fexact", path)
    # Single-phase init registers the module in sys.modules: undo that.
    previous = sys.modules.get("fexact")
    try:
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    finally:
        if previous is None:
            sys.modules.pop("fexact", None)
        else:
            sys.modules["fexact"] = previous
    return module
