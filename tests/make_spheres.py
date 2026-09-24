import math
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox, BRepPrimAPI_MakeCylinder, BRepPrimAPI_MakeSphere
from OCP.BRepAlgoAPI import BRepAlgoAPI_Cut, BRepAlgoAPI_Fuse
from OCP.BRepFilletAPI import BRepFilletAPI_MakeFillet
from OCP.BRepMesh import BRepMesh_IncrementalMesh
from OCP.StlAPI import StlAPI_Writer
from OCP.gp import gp_Ax2, gp_Pnt, gp_Dir
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_EDGE
from OCP.BRepAdaptor import BRepAdaptor_Curve
from OCP.TopoDS import TopoDS
from OCP.GProp import GProp_GProps
from OCP.BRepGProp import BRepGProp
from OCP.BRepTools import BRepTools
def cyl(x,y,z,r,h,d=(0,0,1)): return BRepPrimAPI_MakeCylinder(gp_Ax2(gp_Pnt(x,y,z),gp_Dir(*d)),r,h).Shape()
def sph(x,y,z,r): return BRepPrimAPI_MakeSphere(gp_Pnt(x,y,z),r).Shape()
cut=lambda a,b: BRepAlgoAPI_Cut(a,b).Shape(); fuse=lambda a,b: BRepAlgoAPI_Fuse(a,b).Shape()
s=BRepPrimAPI_MakeBox(70,40,12).Shape()
# fillet all edges at one corner region: fillet all 12 edges r=2 -> spherical corner octants (EXPECT skipped)
f=BRepFilletAPI_MakeFillet(s); ex=TopExp_Explorer(s,TopAbs_EDGE)
while ex.More(): f.Add(2.0,TopoDS.Edge(ex.Current())); ex.Next()
s=f.Shape()
s=fuse(s,cut(sph(12,12,12,5),BRepPrimAPI_MakeBox(gp_Pnt(0,0,0),30,30,12).Shape()))   # hemisphere dome on top
s=cut(s,sph(30,12,12+2.5,4))            # shallow dimple: sphere R4 center 2.5 above top -> depth 1.5
s=fuse(s,cyl(45,12,12,2,6)); s=fuse(s,sph(45,12,18,2))  # capsule pin with hemispherical tip
s=cut(s,cyl(58,12,5,1.5,8)); s=cut(s,sph(58,12,5,1.5))    # ball-ended blind hole
s=cut(s,sph(20,30,6,2.5))              # internal spherical void (closed cavity)
ball=sph(90,20,6,6)                      # separate loose ball
pr=GProp_GProps(); BRepGProp.VolumeProperties_s(s,pr); v1=pr.Mass(); BRepGProp.VolumeProperties_s(ball,pr); print("true volume part",v1,"ball",pr.Mass())
from OCP.TopoDS import TopoDS_Compound
from OCP.BRep import BRep_Builder
comp=TopoDS_Compound(); b=BRep_Builder(); b.MakeCompound(comp); b.Add(comp,s); b.Add(comp,ball)
BRepMesh_IncrementalMesh(comp,0.01,False,0.25,True); w=StlAPI_Writer(); w.ASCIIMode=False; w.Write(comp,"test_sphere.stl")
