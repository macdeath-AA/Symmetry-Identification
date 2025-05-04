from OCC.Core.STEPControl import STEPControl_Reader
from OCC.Display.SimpleGui import init_display
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopAbs import TopAbs_FACE
from OCC.Core.BRepGProp import brepgprop_SurfaceProperties
from OCC.Core.GProp import GProp_GProps
from OCC.Core.BRep import BRep_Tool
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCC.Core.BRepTools import breptools_UVBounds
from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Pln, gp_Dir
from OCC.Core.Geom import Geom_Plane
import numpy as np
import random

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
            'face': face,
            'center': (center.X(), center.Y(), center.Z()),
            'area': area,
            'normal': (normal_vec.X(), normal_vec.Y(), normal_vec.Z())
        })
        explorer.Next()
    
    return face_properties

def vector_between(p1,p2):
    return (p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2])

def normalize(vec):
    norm = np.linalg.norm(vec)
    if norm == 0:
        return vec
    return (vec[0] / norm, vec[1] / norm, vec[2] / norm)

def midpoint(p1, p2):
    return ((p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2, (p1[2] + p2[2]) / 2)

def candidate_planes(face_properties):
    candidate_planes = []

    for i in range(len(face_properties)):
        for j in range(i+1, len(face_properties)):
            f1 = face_properties[i]
            f2 = face_properties[j]

            center1 = f1['center']
            center2 = f2['center']
            mid = midpoint(center1, center2)
            nvec = normalize(vector_between(center1, center2))

            plane = {
                'face1_index': i,
                'face2_index': j,
                'midpoint': mid,
                'plane_normal': nvec
            }

            candidate_planes.append(plane)
    
    return candidate_planes

def visualize_candidate_planes(display, candidate_planes,count=1):
    
    for i, plane_info in enumerate(candidate_planes[:count]):
        origin = gp_Pnt(*plane_info['midpoint'])
        normal = gp_Dir(*plane_info['plane_normal'])
        plane_geom = gp_Pln(origin, normal)
        size = 100
        plane_face = BRepBuilderAPI_MakeFace(plane_geom, -size, size, -size, size).Face()
        display.DisplayShape(plane_face, update=True, color='BLUE1', transparency=0.5)

def visualize_random_candidate_planes(display, candidate_planes, count=5):
    random_planes = random.sample(candidate_planes,count)
    visualize_candidate_planes(display, random_planes, count= count)


if __name__ == "__main__":
    filename= "JAWSLIDING.step"
    shape = read_step_file(filename)

    # Get face properties
    face_data = get_face_properties(shape)
    print("Total Faces Found:", len(face_data))    

    #candidate planes
    candidate_planes_list = candidate_planes(face_data)
    print("Total Candidate Planes Found:", len(candidate_planes_list))
    
    #Visualization
    display, start_display, _, _ = init_display()
    display.DisplayShape(shape, update=True)

    start_display()
