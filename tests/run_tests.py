"""
Regression tests: build test parts in CAD, mesh them to STL, convert them back, and
check the result against the original design.

  python tests/run_tests.py

Each part must come back as a valid solid, within a small volume error of the true
design, with no more faces than expected (fewer faces = more of it rebuilt as true
surfaces).
"""
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from stl2curves import stl_to_solid, count, volume  # noqa: E402
from OCP.BRepCheck import BRepCheck_Analyzer  # noqa: E402
from OCP.TopAbs import TopAbs_FACE  # noqa: E402

# part: (generator, true volume from the CAD model, allowed volume error, max faces)
CASES = {
    "test_fine.stl": ("make_plate.py", 34615.8107, 0.01, 29),     # holes, pins, rounded corners
    "test_coarse.stl": ("make_plate.py", 34615.8107, 0.01, 29),   # same part, coarse mesh
    "test_hard.stl": ("make_hard.py", 35134.2571, 0.05, 31),      # chamfers, slot, crossing / sloped holes
    "test_lid.stl": ("make_lid.py", 57600.4078, 0.05, 45),        # fillet network, tori, cones
    "test_sphere.stl": ("make_spheres.py", 34303.8040, 0.01, 34),  # domes, dimples, balls, corners
}


def main():
    parts = HERE / "parts"
    parts.mkdir(exist_ok=True)
    for script in sorted({g for g, *_ in CASES.values()}):
        subprocess.run([sys.executable, str(HERE / script)], cwd=parts, check=True,
                       stdout=subprocess.DEVNULL)
    failures = 0
    for name, (_, true_vol, allowed, max_faces) in CASES.items():
        t = time.time()
        shape, info = stl_to_solid(parts / name, 0.01, fuse=False)
        vol, faces = volume(shape), count(shape, TopAbs_FACE)
        valid = BRepCheck_Analyzer(shape).IsValid()
        ok = valid and abs(vol - true_vol) <= allowed and faces <= max_faces
        failures += not ok
        print(f"{'PASS' if ok else 'FAIL'}  {name:<16} volume error {vol - true_vol:+.4f} mm^3, "
              f"{faces} faces (max {max_faces}), valid {valid}, {time.time() - t:.1f}s")
    print("all passed" if not failures else f"{failures} failed")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
