import math
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakeCone
from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
from OCP.BRepFilletAPI import BRepFilletAPI_MakeFillet
from OCP.BRepMesh import BRepMesh_IncrementalMesh
from OCP.StlAPI import StlAPI_Writer
from OCP.gp import gp_Ax2, gp_Pnt, gp_Dir, gp_Trsf, gp_Ax1
from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_EDGE
from OCP.BRepAdaptor import BRepAdaptor_Curve
from OCP.GeomAbs import GeomAbs_Circle
from OCP.TopoDS import TopoDS
from OCP.GProp import GProp_GProps
from OCP.BRepGProp import BRepGProp
def cyl(x,y,z,r,h,d=(0,0,1)): return BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(x,y,z),gp_Dir(*d)),r,h).Shape()
def cone(x,y,z,r1,r2,h): return BRepPrimAPI_MakeCone(gp_Ax2(gp_Pnt(x,y,z),gp_Dir(0,0,1)),r1,r2,h).Shape()
cut=lambda a,b: BRepAlgoAPI_Cut(a,b).Shape(); fuse=lambda a,b: BRepAlgoAPI_Fuse(a,b).Shape()
s=BRepPrimAPI_MakeBox(80,40,12).Shape()
# A: hole with 0.6mm chamfer both ends (EXPECT: rebuilt)
s=cut(s,cyl(10,10,-1,2,14)); s=cut(s,cone(10,10,-0.001,2.6,2,0.601)); s=cut(s,cone(10,10,11.4,2,2.6,0.601))
# B: crossing holes (EXPECT: both skipped - ends not flat)
s=cut(s,cyl(30,-1,6,1.5,42,(0,1,0))); s=cut(s,cyl(25,20,6,1.5,10,(1,0,0)))
# C: slot / obround (EXPECT: skipped)
s=cut(s,cyl(50,10,-1,2,14)); s=cut(s,cyl(56,10,-1,2,14)); s=cut(s,BRepPrimAPI_MakeBox(gp_Pnt(50,8,-1),6,4,14).Shape())
# D: pin with chamfered tip (EXPECT: rebuilt, tip cone stays faceted)
s=fuse(s,cyl(15,30,12,2.5,5)); s=fuse(s,cone(15,30,17,2.5,1.9,0.6))
# F: pin with filleted base
p=cyl(35,30,11,2,7); s=fuse(s,p)
f=BRepFilletAPI_MakeFillet(s); ex=TopExp_Explorer(s,TopAbs_EDGE)
while ex.More():
    e=TopoDS.Edge(ex.Current()); c=BRepAdaptor_Curve(e)
    if c.GetType()==GeomAbs_Circle and abs(c.Circle().Radius()-2)<1e-6 and abs(c.Circle().Location().Z()-12)<1e-6 and abs(c.Circle().Location().X()-35)<1e-6: f.Add(0.8,e)
    ex.Next()
s=f.Shape()
# G: hole through sloped face: cut a wedge off corner then vertical hole there (EXPECT skipped)
w=BRepPrimAPI_MakeBox(gp_Pnt(60,20,0),30,30,30).Shape()
t=gp_Trsf(); t.SetRotation(gp_Ax1(gp_Pnt(60,20,12),gp_Dir(0,1,0)),math.radians(-30)); w=BRepBuilderAPI_Transform(w,t).Shape()
s=cut(s,w); s=cut(s,cyl(66,30,-1,1.8,20))
# H: plain hole deep blind horizontal from side (EXPECT rebuilt)
s=cut(s,cyl(-1,34,5,1.25,20,(1,0,0)))
pr=GProp_GProps(); BRepGProp.VolumeProperties_s(s,pr); print("true volume",pr.Mass())
BRepMesh_IncrementalMesh(s,0.01,False,0.3,True); w=StlAPI_Writer(); w.ASCIIMode=False; w.Write(s,"test_hard.stl")
