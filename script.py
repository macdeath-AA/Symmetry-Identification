from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Display.SimpleGui import init_display
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE
from OCC.Core.BRepGProp import brepgprop_SurfaceProperties
from OCC.Core.GProp import GProp_GProps
from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepTools import breptools_UVBounds
from OCC.Core.gp import gp_Pnt, gp_Vec
import numpy as np

def read_step_file(filename):
    reader = STEPControl_Reader()
    reader.ReadFile(filename)
    reader.TransferRoots()
    return reader.OneShape()

def get_face_properties(shape):
    explorer = TopExp_Explorer(shape, TopAbs_FACE)
    face_properties = []

    while explorer.More():
        face = explorer.Current()
        props = GProp_GProps()
        brepgprop_SurfaceProperties(face, props)
        center = props.CentreOfMass()
        area = props.Mass()
        umin, umax, vmin, vmax = breptools_UVBounds(face)
        normal_vec = gp_Vec()
        BRep_Tool.Surface(face).D1((umin+ umax)/2, (vmin+vmax)/2, gp_Pnt(), normal_vec, gp_Vec())
        face_properties.append({
            'center': (center.X(), center.Y(), center.Z()),
            'area': area,
            'normal': (normal_vec.X(), normal_vec.Y(), normal_vec.Z())
        })
        explorer.Next()
    
    return face_properties

if __name__ == "__main__":
    filename= "GearShaft.step"
    shape = read_step_file(filename)

    # Get face properties
    face_data = get_face_properties(shape)
    print("Total Faces Found:", len(face_data))

    for i, face in enumerate(face_data):
        print(f"Face {i+1}: Center = {face['center']}, Area = {face['area']:.4f}, Normal = {face['normal']}")

    #Visualization
    display, start_display, _,_ = init_display()
    display.DisplayShape(shape, update=True)
    start_display()
