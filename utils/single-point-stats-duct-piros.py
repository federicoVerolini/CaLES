import numpy as np

filename = 'plotyz_Retau1000.dat'
data_avg = np.loadtxt(filename, skiprows=20)

n1 = int(np.sqrt(data_avg.shape[0]))
n2 = int(np.sqrt(data_avg.shape[0]))
data_avg  = np.reshape(data_avg,(n2,n1,14),order='C')
#
# bisector statistics
#
zc_cl = data_avg[n2-1,:,0]+1.0
yc_cl = data_avg[n2-1,:,1]+1.0
fname = 'stats-single-point-duct-centerline.out'
with open(fname, 'w') as file:
    np.savetxt(file, np.c_[zc_cl,yc_cl,data_avg[n2-1,:,2:]], \
               fmt='%17.9e', delimiter='')
#
# diagonal statistics, uniform grid in the cross-section
#
zc_diag = np.diag(data_avg[:,:,0])+1.0
yc_diag = np.diag(data_avg[:,:,1])+1.0
avg_diag = np.zeros((n2,14))
for k in range(2,14):
    avg_diag[:,k] = np.diag(data_avg[:,:,k])
fname = 'stats-single-point-duct-diagonal.out'
with open(fname, 'w') as file:
    np.savetxt(file, np.c_[zc_diag,yc_diag,avg_diag[:,2:]], \
               fmt='%16.6e', delimiter='')