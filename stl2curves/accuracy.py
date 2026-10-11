"""How far a rebuilt surface may stray from the mesh: one setting for the trade between a
faithful part and a finished one (a STEP whose every curve is a surface but half a
millimetre off can be worth more than one accurate to a micron and half triangles).

--accuracy MM (or the environment variable STL2CURVES_ACCURACY, which worker processes
inherit) scales every limit on the surfaces fitted where no exact design surface is sure:
smooth blends and freeform surfaces (how close to the corners, how far they may bow
between them), fillets on rough meshes, and the coarse-export cylinders and cones. Exact
surfaces on clean meshes are as true as the mesh whatever it is. The limits are written
as tuned at TUNED (0.02 mm); the default, 0.05 mm (2026-10-11), scales them by 2.5: a mesh
never holds an exact radius and a CAD program's own spline export strays a few
hundredths, so held to 0.02 real shapes were left as triangles and patches (quick set
3,846 -> 3,352 faces, the measured deviation from the mesh no worse: 99% of points
within 0.02 mm).
"""
import os
import sys

TUNED = 0.02        # mm: the accuracy every limit was tuned at (as written in the modules)
DEFAULT = 0.05      # mm: the accuracy used unless one is asked for
# the limits it scales, as (module, name); some are copies of others' (imported by name)
LIMITS = [("blends", "MAX_DEVIATION"), ("blends", "MAX_EDGE_GAP"),
          ("build", "MAX_DEVIATION"), ("build", "MAX_EDGE_GAP"), ("build", "HUG_SLACK"),
          ("freeform", "FIT_DEV"), ("freeform", "BULGE_MM"), ("pipes", "FIT_DEV"), ("extrude", "FIT_DEV"),
          ("features", "NOISY_MAX"), ("features", "NOISY_TYPICAL"),
          ("fillets", "FILLET_TOL"), ("fillets", "GROW_TOL"),
          ("threads", "PIECE_TOL"), ("features", "CORNER_FIT_DEV")]


def mm():
    """The accuracy asked for (mm)."""
    return float(os.environ.get("STL2CURVES_ACCURACY") or DEFAULT)


SCALE = mm() / TUNED        # (the modules' limits are their tuned values times this)


def set_accuracy(value):
    """Use this accuracy (mm) from now on, here and in worker processes started after
    (stop running ones first: workers.finish)."""
    global SCALE
    value = float(value)
    if value <= 0:
        raise ValueError("accuracy must be more than 0 mm")
    os.environ["STL2CURVES_ACCURACY"] = repr(value)
    new = value / TUNED
    for module, name in LIMITS:
        m = sys.modules.get(f"{__package__}.{module}")
        if m is not None and hasattr(m, name):
            setattr(m, name, getattr(m, name) / SCALE * new)
    SCALE = new
