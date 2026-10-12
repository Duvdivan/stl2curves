"""Transition fillets: where rounded edges of different radii meet, as CAD makes them.

test_corners.stl: a box with its upright edges rounded r 3, then its top edge loop r 1
(the ball rolls on round the r 3 arcs: cylinders and tori) and its bottom loop r 4 (a
ball bigger than the arcs can't roll round them: the corners are spline patches).
test_corners2.stl: a box with all its edges rounded in one go, top r 1, upright r 2,
bottom r 3 (each corner a spline vertex blend between three radii)."""
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.BRepFilletAPI import BRepFilletAPI_MakeFillet
from OCP.BRepMesh import BRepMesh_IncrementalMesh
from OCP.StlAPI import StlAPI_Writer
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_EDGE
from OCP.TopTools import TopTools_ShapeMapHasher
from OCP.BRepAdaptor import BRepAdaptor_Curve
from OCP.GeomAbs import GeomAbs_Line
from OCP.TopoDS import TopoDS
from OCP.GProp import GProp_GProps
from OCP.BRepGProp import BRepGProp


def edges(s, pred):
    out, seen = [], set()
    ex = TopExp_Explorer(s, TopAbs_EDGE)
    while ex.More():
        e = TopoDS.Edge(ex.Current())
        ex.Next()
        key = TopTools_ShapeMapHasher()(e)
        if key in seen:
            continue
        seen.add(key)
        c = BRepAdaptor_Curve(e)
        p0, p1 = c.Value(c.FirstParameter()), c.Value(c.LastParameter())
        if pred(c, p0, p1):
            out.append(e)
    return out


def fillet(s, items):
    f = BRepFilletAPI_MakeFillet(s)
    for r, e in items:
        f.Add(r, e)
    return f.Shape()


def write(s, name):
    pr = GProp_GProps()
    BRepGProp.VolumeProperties_s(s, pr)
    print(name, "true volume", pr.Mass())
    BRepMesh_IncrementalMesh(s, 0.01, False, 0.25, True)
    w = StlAPI_Writer()
    w.ASCIIMode = False
    w.Write(s, name)


H = 14
upright = lambda c, p0, p1: c.GetType() == GeomAbs_Line and abs(p0.Z() - p1.Z()) > 1
at = lambda z: (lambda c, p0, p1: abs(p0.Z() - z) < 1e-6 and abs(p1.Z() - z) < 1e-6)

s = BRepPrimAPI_MakeBox(60, 40, H).Shape()
s = fillet(s, [(3.0, e) for e in edges(s, upright)])
s = fillet(s, [(1.0, e) for e in edges(s, at(H))])
s = fillet(s, [(4.0, e) for e in edges(s, at(0))])
write(s, "test_corners.stl")

s = BRepPrimAPI_MakeBox(40, 30, 12).Shape()
s = fillet(s, [(1.0, e) for e in edges(s, at(12))] + [(2.0, e) for e in edges(s, upright)]
           + [(3.0, e) for e in edges(s, at(0))])
write(s, "test_corners2.stl")
