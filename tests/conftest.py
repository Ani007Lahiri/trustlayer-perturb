"""Test-suite config. RP01 belt-and-suspenders: pin PYTHONHASHSEED so any residual
built-in hash() use is reproducible. The primary fix is _stable_hash() in conformal.py
and trust.py; this guarantees determinism for the whole test run regardless."""
import os
import sys

os.environ.setdefault("PYTHONHASHSEED", "0")

# Re-exec once with the pinned hash seed if it wasn't set before interpreter start
# (PYTHONHASHSEED only takes effect at process launch).
if os.environ.get("_TRUSTLAYER_HASHSEED_REEXEC") != "1" and os.environ.get("PYTHONHASHSEED") != "0":
    os.environ["_TRUSTLAYER_HASHSEED_REEXEC"] = "1"
    os.environ["PYTHONHASHSEED"] = "0"
    os.execv(sys.executable, [sys.executable] + sys.argv)
