"""How far a rebuilt surface may stray from the mesh: one setting for the trade between a
faithful part and a finished one (a STEP whose every curve is a surface but half a
millimetre off can be worth more than one accurate to a micron and half triangles).

--accuracy MM (or the environment variable STL2CURVES_ACCURACY, which worker processes
inherit) scales every limit on the surfaces fitted where no exact design surface is sure:
smooth blends and freeform surfaces (how close to the corners, how far they may bow
between them), fillets on rough meshes, and the coarse-export cylinders and cones. Exact
surfaces on clean meshes are as true as the mesh whatever it is. The default, 0.02 mm,
leaves every limit as tuned.
"""
import os
import sys

DEFAULT = 0.02      # mm: the accuracy every limit was tuned at
# the limits it scales, as (module, name); some are copies of others' (imported by name)
LIMITS = [("blends", "MAX_DEVIATION"), ("blends", "MAX_EDGE_GAP"),
          ("build", "MAX_DEVIATION"), ("build", "MAX_EDGE_GAP"), ("build", "HUG_SLACK"),
          ("freeform", "FIT_DEV"), ("freeform", "BULGE_MM"), ("pipes", "FIT_DEV"), ("extrude", "FIT_DEV"),
          ("features", "NOISY_MAX"), ("features", "NOISY_TYPICAL"),
          ("fillets", "FILLET_TOL"), ("fillets", "GROW_TOL")]


def mm():
    """The accuracy asked for (mm)."""
    return float(os.environ.get("STL2CURVES_ACCURACY") or DEFAULT)


SCALE = mm() / DEFAULT      # (the modules' limits are their tuned values times this)


def set_accuracy(value):
    """Use this accuracy (mm) from now on, here and in worker processes started after
    (stop running ones first: workers.finish)."""
    global SCALE
    value = float(value)
    if value <= 0:
        raise ValueError("accuracy must be more than 0 mm")
    os.environ["STL2CURVES_ACCURACY"] = repr(value)
    new = value / DEFAULT
    for module, name in LIMITS:
        m = sys.modules.get(f"{__package__}.{module}")
        if m is not None and hasattr(m, name):
            setattr(m, name, getattr(m, name) / SCALE * new)
    SCALE = new
