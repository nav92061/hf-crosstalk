"""hf_crosstalk — importable access to the six analysis skills.

The skills live in ``skills/<name>/kernel.py`` with hyphenated directory names,
which are not valid Python module names. This shim loads them by file path and
exposes each under an underscored alias::

    from hf_crosstalk import crosstalk_score as ck
    floor = ck.housekeeping_floor(expr, "GAPDH", 0.005)

or dynamically::

    import hf_crosstalk
    ck = hf_crosstalk.load("crosstalk-score")
    print(hf_crosstalk.SKILLS)

Every kernel is pure stdlib + numpy/pandas/scipy (plus ``requests`` for
depmap-lineage) and references no agent-platform API, so this package works
anywhere Python does.
"""

import importlib.util as _ilu
import os as _os
import sys as _sys

__version__ = "1.0.0"

SKILLS = ("geo-bulk-de", "hpa-secretome", "cellchat-lr",
          "tcga-pancan", "crosstalk-score", "depmap-lineage")

_HERE = _os.path.dirname(_os.path.abspath(__file__))
_CANDIDATES = (_os.path.join(_HERE, "skills"),
               _os.path.join(_os.path.dirname(_HERE), "skills"))
_CACHE = {}


def skills_dir():
    """Return the directory holding the skill subdirectories.

    Works both from a source checkout and from an installed wheel, where the
    skills are shipped as package data.
    """
    for path in _CANDIDATES:
        if _os.path.isdir(path):
            return path
    raise FileNotFoundError(
        "skills/ not found; looked in: %s" % ", ".join(_CANDIDATES))


def skill_path(name):
    """Absolute path to a skill's directory. Raises ValueError if unknown."""
    if name not in SKILLS:
        raise ValueError("unknown skill %r; choose from %s"
                         % (name, ", ".join(SKILLS)))
    return _os.path.join(skills_dir(), name)


def read_doc(name):
    """Return the SKILL.md text for ``name``."""
    with open(_os.path.join(skill_path(name), "SKILL.md"), encoding="utf-8") as fh:
        return fh.read()


def load(name):
    """Import and return a skill's kernel module, caching the result."""
    if name in _CACHE:
        return _CACHE[name]
    kernel = _os.path.join(skill_path(name), "kernel.py")
    if not _os.path.exists(kernel):
        raise FileNotFoundError("no kernel.py for skill %r at %s" % (name, kernel))
    mod_name = "hf_crosstalk._k_" + name.replace("-", "_")
    spec = _ilu.spec_from_file_location(mod_name, kernel)
    module = _ilu.module_from_spec(spec)
    _sys.modules[mod_name] = module
    spec.loader.exec_module(module)
    _CACHE[name] = module
    return module


def __getattr__(attr):
    """Expose underscored skill aliases lazily (PEP 562).

    ``hf_crosstalk.crosstalk_score`` loads ``skills/crosstalk-score/kernel.py``
    only when first accessed, so importing this package stays cheap.
    """
    hyphenated = attr.replace("_", "-")
    if hyphenated in SKILLS:
        return load(hyphenated)
    raise AttributeError("module %r has no attribute %r" % (__name__, attr))


def __dir__():
    return sorted(list(globals()) + [s.replace("-", "_") for s in SKILLS])
