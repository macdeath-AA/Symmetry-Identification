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
from scipy.spatial import cKDTree
import numpy as np
import random

def readStepFile(filename):
    reader = STEPControl_Reader()
    reader.ReadFile(filename)
    reader.TransferRoots()
    return reader.OneShape()

def getFaceProperties(shape):
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

def vectorBetween(p1,p2):
    return (p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2])

def normalize(vec):
    norm = np.linalg.norm(vec)
    if norm == 0:
        return vec
    return (vec[0] / norm, vec[1] / norm, vec[2] / norm)

def midpoint(p1, p2):
    return ((p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2, (p1[2] + p2[2]) / 2)

def candidatePlanes(face_properties, area_tol=0.9, normal_tol=0.8):
    candidate_planes = []
    n = len(face_properties)

    for i in range(n):
        ci, ni, Ai = face_properties[i]['center'], face_properties[i]['normal'], face_properties[i]['area']
        
        for j in range(i+1,n):
            cj,nj, Aj = face_properties[j]['center'], face_properties[j]['normal'], face_properties[j]['area']
            
            #1. area ratio check
            if min(Ai, Aj) / max(Ai, Aj) < area_tol:
                continue
            
            #2. normal alignment check
            if abs(np.dot(ni, nj )) < normal_tol:
                continue

            nvec = normalize(vectorBetween(ci, cj))            

            if not bisectNormals(ni, nj, nvec):
                continue

            mid = midpoint(ci, cj)            

            plane = {
                'face1_index': i,
                'face2_index': j,
                'midpoint': mid,
                'plane_normal': nvec
            }

            candidate_planes.append(plane)
    
    return candidate_planes

def visualizeCandidatePlanes(display, candidate_planes,count=1):
    
    for i, plane_info in enumerate(candidate_planes[:count]):
        origin = gp_Pnt(*plane_info['midpoint'])
        normal = gp_Dir(*plane_info['plane_normal'])
        plane_geom = gp_Pln(origin, normal)
        size = 100
        plane_face = BRepBuilderAPI_MakeFace(plane_geom, -size, size, -size, size).Face()
        display.DisplayShape(plane_face, update=True, color='BLUE1', transparency=0.5)

def reflectPoint(point, plane_point, plane_normal):
    v = np.array(point) - plane_point
    return plane_point + v -2 * np.dot(v,plane_normal) * plane_normal

def reflectVector(vec, plane_normal):
    return vec - 2 * np.dot(vec, plane_normal) * plane_normal

def bisectNormals(n1,n2, plane_normal, tol = 1e-2):
    n1 = np.array(n1)
    n2 = np.array(n2)
    return np.abs(np.dot(n1,plane_normal) + np.dot(n2,plane_normal)) < tol

def PairwiseMatchingScore(face_data, plane):
    pts = np.array([f['center'] for f in face_data])
    norms = np.array([f['normal'] for f in face_data])
    areas = np.array([f['area'] for f in face_data])
    tree = cKDTree(pts)

    p0 = np.array(plane['midpoint'])
    N = np.array(plane['plane_normal'])

    dist_error = 0.0
    norm_error = 0.0
    area_error = 0.0
    M = len(face_data)

    for i in range(M):
        c,n,a = pts[i], norms[i], areas[i]

        #reflect center --> find nearest real face
        c_ref = reflectPoint(c, p0, N)
        d,j = tree.query(c_ref)
        dist_error += d 

        #reflect normal --> compare alignment
        n_ref = reflectVector(n, N)
        n_ref /= np.linalg.norm(n_ref) 
        n_j = norms[j]/ np.linalg.norm(norms[j])
        norm_error += (1 - abs(np.dot(n_ref, n_j)))

        #compare areas
        area_error += abs(a - areas[j]) / max(a, areas[j])
    
    return {
        'avg_dist_error': dist_error / M,
        'avg_norm_error': norm_error / M,
        'avg_area_error': area_error / M
    }

def findBestPlane(face_data, candidate_planes):
    best, best_score = None, None

    for plane in candidate_planes:
        errs = PairwiseMatchingScore(face_data, plane)

        total_err = errs['avg_dist_error'] + errs['avg_norm_error'] + errs['avg_area_error']

        if best is None or total_err < best_score:
            best_score, best = total_err, plane
    
    return best


if __name__ == "__main__":
    filename= "JAWSLIDING.step"
    shape = readStepFile(filename)

    # Get face properties
    face_data = getFaceProperties(shape)
    print("Total Faces Found:", len(face_data))    

    #candidate planes
    candidate_planes_list = candidatePlanes(face_data)
    print("Total Candidate Planes Found:", len(candidate_planes_list))

    #best plane
    best_plane = findBestPlane(face_data, candidate_planes_list)
    print("Best Plane Found:", best_plane)
    
    #Visualization
    display, start_display, _, _ = init_display()
    display.DisplayShape(shape, update=True)

    visualizeCandidatePlanes(display, [best_plane], count=11)

    start_display()
