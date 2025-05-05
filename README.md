# Symmetry Identification with Performance Constraints
This project detects significant symmetry in 3D models (specifically STEP files) by analyzing the faces of the model. The tool works by identifying potential mirror planes that could reflect the shape symmetrically. If no significant symmetry is found, it provides a message indicating that.

## Installation
Install the required Python libraries: OCC, numpy, SciPy, FreeCAD (optional)
To install [pyhtonocc](https://github.com/tpaviot/pythonocc-core) follow these steps.

   ```
   conda create --name=pyoccenv python=3.9
conda activate pyoccenv
conda install -c conda-forge pythonocc-core=7.9.0
   ```

To install the required Python libraries, run the following command in your environment. 

<pre>
   ```pip install numpy scipy```
</pre>

## Usage
You can run the script directly from the command line using Python.
<pre>
   ``` python script.py```
</pre>

Make sure the STEP files are in the same directory as ``script.py``

## Methodology Overview
1. **Face Extraction**: All faces of the shape are processed to compute their centroids, areas, and normal vectors.
2. **Candidate Plane Generation**: A set of potential symmetry planes is constructed, typically using face-based heuristics (e.g. average planes between similar normals).
3. **Scoring Planes**: Each candidate plane is evaluated for symmetry based on how well it aligns matching face pairs across the plane. The score considers:
   - Centroid Distance Error
   - Normal Vector Alignment Error
   - Area Error
4. **Selecting Distinct Planes**: Top planes with the lowest errors are selected. To ensure uniqueness, only planes with non-parallel normals are included (based on dot product threshold).
5. **Visualization**: The original shape and the selected symmetry planes are rendered in a 3D window.

   

