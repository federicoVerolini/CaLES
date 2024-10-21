import math
import os
import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk

# generating the .vtk file of the Y-Z averaged section
data=np.loadtxt("velMean_YZ.out")
coords=data[:,:2]
fields=data[:,2:]
Ny=int(math.sqrt(coords.shape[0]))
Nz=int(math.sqrt(coords.shape[0]))
Nx=1

points=vtk.vtkPoints()
for coord in coords:
    points.InsertNextPoint(coord[0],coord[1],0.0)
structured_grid=vtk.vtkStructuredGrid()
structured_grid.SetDimensions(Ny,Nz,Nx)
structured_grid.SetPoints(points)

campi=["um","vm","wm","um2","vm2","wm2"]
for i in range(fields.shape[1]):
    field_data=numpy_to_vtk(fields[:,i],deep=True)
    field_data.SetName(f"{campi[i]}")
    structured_grid.GetPointData().AddArray(field_data)

writer=vtk.vtkStructuredGridWriter()
writer.SetFileName("velMean_YZ.vtk")
writer.SetInputData(structured_grid)
writer.Write()