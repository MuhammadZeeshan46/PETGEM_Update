import meshio
import numpy as np
from pathlib import Path
from petsc4py import PETSc

# ------------------------------------------------------------------------------
#              USER PARAMS
# ------------------------------------------------------------------------------
mesh_file = "mesh.msh"
sigma_x = np.array([1., 0.01, 1., 3.3333], dtype=float)
sigma_y = np.array([1., 0.01, 1., 3.3333], dtype=float)
sigma_z = np.array([1., 0.01, 1., 3.3333], dtype=float)

# ------------------------------------------------------------------------------
#              IMPORT MESH
# ------------------------------------------------------------------------------
mesh = meshio.read(mesh_file)

# Determine number of materials
num_materials = len(mesh.cells)

# Number of elements, dimensions and connectivity
nElems = 0
for i in np.arange(num_materials):
	elemsN = mesh.cells[i].data
	nElems += np.shape(elemsN)[0]

# Number of nodes and coordinates
nodes = mesh.points
nnodes = np.shape(nodes)[0]

# Print mesh info
print("Number of elements:", nElems)
print("Number of nodes:", nnodes)
print("Mesh points:")
print("\t Min x-coordinate:", np.min(nodes[:,0]))
print("\t Max x-coordinate:", np.max(nodes[:,0]))
print("\t Min y-coordinate:", np.min(nodes[:,1]))
print("\t Max y-coordinate:", np.max(nodes[:,1]))
print("\t Min z-coordinate:", np.min(nodes[:,2]))
print("\t Max z-coordinate:", np.max(nodes[:,2]))

# ------------------------------------------------------------------------------
#              COMPUTE RESISTIVITY MODEL
# ------------------------------------------------------------------------------
sigma = np.zeros((nElems,3), dtype=float)
elemsS = np.copy(mesh.cell_data_dict["gmsh:physical"]["tetra"])
elemsS -= 1

for i in np.arange(nElems):
	sigma[i,0] = sigma_x[elemsS[i]]	# Isotropic model
	sigma[i,1] = sigma_y[elemsS[i]] # Isotropic model
	sigma[i,2] = sigma_z[elemsS[i]] # Isotropic model

elemsN = np.zeros((nElems,4), dtype=int)
indx = 0

for i in np.arange(num_materials):
	tmp = mesh.cells[i].data
	nElems_local = np.shape(tmp)[0]
	for j in np.arange(nElems_local):
		elemsN[indx,: ] = tmp[j,:]
		indx += 1

# ------------------------------------------------------------------------------
#  EXPORT MESH 
# ------------------------------------------------------------------------------
mesh_petgem = meshio.Mesh(points=mesh.points,
	                      cells=[("tetra", elemsN)], 
 	                      cell_data={"gmsh:physical":[mesh.cell_data_dict["gmsh:physical"]["tetra"]], 
                                     "gmsh:geometrical":[mesh.cell_data_dict["gmsh:geometrical"]["tetra"]],
                                     "x_conductivity": [sigma[:,0]], 
 	                                 "y_conductivity": [sigma[:,1]], 
 	                                 "z_conductivity": [sigma[:,2]]})

# write mesh
tmp_file_name = Path(mesh_file)
out_mesh_filename = tmp_file_name.stem + '_resistivity.vtu'
mesh_petgem.write(out_mesh_filename, file_format='vtu')

# Export resistivity model for x-component
vector = PETSc.Vec().createWithArray(sigma[:,0], comm=PETSc.COMM_SELF)
vector.setUp()
viewer = PETSc.Viewer().createBinary('resistivity_data_x.dat', mode='w', comm=PETSc.COMM_SELF)
vector.view(viewer)

# Destroy the vector and viewer
vector.destroy()
viewer.destroy()


# Export resistivity model for y-component
vector = PETSc.Vec().createWithArray(sigma[:,1], comm=PETSc.COMM_SELF)
vector.setUp()
viewer = PETSc.Viewer().createBinary('resistivity_data_y.dat', mode='w', comm=PETSc.COMM_SELF)
vector.view(viewer)

# Destroy the vector and viewer
vector.destroy()
viewer.destroy()

# Export resistivity model for z-component
vector = PETSc.Vec().createWithArray(sigma[:,2], comm=PETSc.COMM_SELF)
vector.setUp()
viewer = PETSc.Viewer().createBinary('resistivity_data_z.dat', mode='w', comm=PETSc.COMM_SELF)
vector.view(viewer)

# Destroy the vector and viewer
vector.destroy()
viewer.destroy()

# ------------------------------------------------------------------------------
#              EXPORT RECEIVERS FILE
# ------------------------------------------------------------------------------
# Import receiver coordinates
receivers = np.loadtxt('receivers.txt')
receivers = receivers.flatten()

vector = PETSc.Vec().createWithArray(receivers, comm=PETSc.COMM_SELF)
vector.setUp()
viewer = PETSc.Viewer().createBinary('receivers.dat', mode='w', comm=PETSc.COMM_SELF)
vector.view(viewer)

# Destroy the vector and viewer
vector.destroy()
viewer.destroy()
