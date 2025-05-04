from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Display.SimpleGui import init_display
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE
from OCC.Core.BRepGProp import brepgprop_SurfaceProperties
from OCC.Core.GProp import GProp_GProps

reader = STEPControl_Reader()
status = reader.ReadFile("GearShaft.step")
reader.TransferRoots()
shape = reader.OneShape()

display,start_display,_,_ = init_display()
display.DisplayShape(shape, update=True)

face_explorer = TopExp_Explorer(shape, TopAbs_FACE)
face_centers = []

while face_explorer.More():
    face = face_explorer.Current()
    props = GProp_GProps()
    brepgprop_SurfaceProperties(face, props)
    
    # Get the center of the face
    center = props.CentreOfMass()
    face_centers.append((center.X(), center.Y(), center.Z()))
    face_explorer.Next()

print("Face centers:")
for i, c in enumerate(face_centers):
    print(f"Face {i + 1}: Center = {c}")

start_display()


