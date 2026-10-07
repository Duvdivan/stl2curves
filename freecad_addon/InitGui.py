# FreeCAD runs this file (with exec, so it has no __file__) when the GUI starts; the
# add-on's code lives in stl2curves_freecad.py, which FreeCAD's loader put on sys.path.
import stl2curves_freecad

stl2curves_freecad.setup()
