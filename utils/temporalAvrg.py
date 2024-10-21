import math
import os
import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk

matrices=[]

folder_path='data/'
for filename in os.listdir(folder_path):
    file_path=os.path.join(folder_path,filename)
    matrix=[]
    with open(file_path, 'r') as file:
        for line in file:
            row=line.rstrip('\n')
            row=row.split()
            row=np.array(row,dtype=float) 
            matrix.append(row)
    matrices.append(matrix)

matrices=np.array(matrices)

avrg_matrix=np.zeros(matrices[0].shape,dtype=float)
for k in range(len(matrices)):
    for i in range(matrices[k].shape[0]):
        for j in range(matrices[k].shape[1]):
            avrg_matrix[i][j]=avrg_matrix[i][j]+matrices[k][i][j]
        
for i in range(avrg_matrix.shape[0]):
    for j in range(avrg_matrix.shape[1]):        
        avrg_matrix[i][j]=avrg_matrix[i][j]/len(matrices)

output_file='velMean_YZ.out'
with open(output_file,'w') as file:
    file.write(f"# Y: Z: um: vm: wm: um2: vm2: wm2:\n")
    for i in range(avrg_matrix.shape[0]):
        for j in range(avrg_matrix.shape[1]):
            if j==avrg_matrix.shape[1]-1:
                file.write(f"{avrg_matrix[i][j]:.7e}")
            else:
                file.write(f"{avrg_matrix[i][j]:.7e} ")
        if i!=avrg_matrix.shape[0]-1:
            file.write("\n")
'''
output_file='velMean_YZ.csv'
with open(output_file,'w') as file:
    file.write(f"Y, Z, um, vm, wm, um2, vm2, wm2,\n")
    for i in range(avrg_matrix.shape[0]):
        for j in range(avrg_matrix.shape[1]):
            if j==avrg_matrix.shape[1]-1:
                file.write(f"{avrg_matrix[i][j]:.7e}")
            else:
                file.write(f"{avrg_matrix[i][j]:.7e}, ")
        if i!=avrg_matrix.shape[0]-1:
            file.write("\n")
'''
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