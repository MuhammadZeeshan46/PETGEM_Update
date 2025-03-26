import numpy as np
import matplotlib.pyplot as plt
import scipy
from petsc4py import PETSc


# ------------------------------------------------------------------------------
#              USER PARAMS
# ------------------------------------------------------------------------------
fields_filename = "output/Ex.dat"
receivers_filename = "receivers.dat"

# ------------------------------------------------------------------------------
#              IMPORT FIELDS AND RECEIVERS
# ------------------------------------------------------------------------------
# Create vector 
vec = PETSc.Vec().create(PETSc.COMM_SELF)
    
# Load the vector from a binary file
viewer = PETSc.Viewer().createBinary(fields_filename, mode='r')
vec.load(viewer)

# Convert the PETSc vector to a NumPy array
Ex = vec.getArray()

# Destroy PETSc objects
vec.destroy()
viewer.destroy()

# Load the receiver coordinates 
viewer = PETSc.Viewer().createBinary(receivers_filename, mode='r')
vec.load(viewer)

# Convert the PETSc vector to a NumPy array
receivers = vec.getArray()
# Extract only x-coordinates 
x_coordinates = receivers[::3]


# Load reference
reference = scipy.io.loadmat('reference.mat')
reference = reference['Reference']

matlab = scipy.io.loadmat('results.mat')
matlab = matlab['Ex']




v1 = np.abs(reference)
v2 = np.abs(Ex)

#v3 = np.abs(matlab)
#print(v2)


# Plot Ex component
plt.figure(figsize=(10, 6))
    
plt.semilogy(x_coordinates, v1, 'r-', marker='+')
plt.semilogy(x_coordinates, v2, 'b-', marker='o')

#plt.semilogy(x_coordinates, v3, 'g-', marker='s')
plt.title('Ex field')
plt.xlabel('Offset')
plt.ylabel('|Ex|')
plt.legend(['Reference', 'PETSC version'])

   
plt.tight_layout()
plt.show()
    


    
