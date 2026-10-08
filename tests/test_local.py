"""Check stl2curves.local against OCC's whole-shape ReShape: the same result, quicker.

    python tests/test_local.py [part.stl]     (default: tests/parts/test_hard.stl)

Builds the part from its bare facets, swaps one edge for a copy of itself through
apply_local and through a plain ReShape, and compares faces, open edges and validity.
"""
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

from OCP.BRepBuilderAPI import BRepBuilderAPI_Copy  # noqa: E402
from OCP.BRepCheck import BRepCheck_Analyzer  # noqa: E402
from OCP.BRepTools import BRepTools_ReShape  # noqa: E402
from OCP.TopoDS import TopoDS  # noqa: E402

from stl2curves import build, local  # noqa: E402
from stl2curves.features import analyze, load_stl  # noqa: E402


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else str(HERE / "parts" / "test_hard.stl")
    pts, tris = load_stl(path)[:2]
    mesh = analyze(pts, tris)[0]
    comp, shells, _ = build.build_faces(mesh, [], 1e-6)
    shape, _ = build.sew(comp, 0.01, shells)
    m = local.edge_faces(shape)
    print(f"{path}: {len(local.faces_of(shape))} faces, {m.Extent()} edges")
    ok = True
    for i in (1, m.Extent() // 2, m.Extent()):
        edge = TopoDS.Edge(m.FindKey(i))
        twin = TopoDS.Edge(BRepBuilderAPI_Copy(edge).Shape())
        reshape = BRepTools_ReShape()
        reshape.Replace(edge, twin)
        t = time.perf_counter()
        whole = reshape.Apply(shape)
        t_whole = time.perf_counter() - t
        reshape = BRepTools_ReShape()
        reshape.Replace(edge, twin)
        t = time.perf_counter()
        part, _ = local.apply_local(shape, reshape, [TopoDS.Face(f) for f in m.FindFromIndex(i)])
        t_local = time.perf_counter() - t
        a, b = local.edge_faces(whole), local.edge_faces(part)
        same = (len(local.faces_of(whole)) == len(local.faces_of(part)) and a.Extent() == b.Extent()
                and len(local.open_edges(a)) == len(local.open_edges(b)) == 0
                and BRepCheck_Analyzer(part).IsValid() == BRepCheck_Analyzer(whole).IsValid()
                and b.Contains(twin) and not b.Contains(edge))
        ok &= same
        print(f"  edge {i}: whole {1000 * t_whole:.1f} ms, local {1000 * t_local:.1f} ms, same result {same}")
    print("all passed" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
