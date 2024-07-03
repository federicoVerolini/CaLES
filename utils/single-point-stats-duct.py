from scipy.interpolate import interp1d
import numpy as np
import sys
import os
import re
import glob
#
# compute mean pressure gradient
#
use_new_arrays = False
shuffle_data   = False
args = sys.argv[1:]
casedir = args[3] if len(args) > 3 else ''
datadir = casedir+''
print('Data directory: ' + datadir+'\n')
resultsdir = casedir+'results/'
os.makedirs(resultsdir,exist_ok=True)
#
# import parameters
#
sys.path.append(casedir) # put input.py in the case directory
from input import *
visc = visci**(-1)
reb = 2*h*ub/visc
tbeg   = float(args[0]) if len(args) > 0 else tbeg
tend   = float(args[1]) if len(args) > 1 else tend
fldstp = int(  args[2]) if len(args) > 2 else fldstp
fname_dpdx = datadir+'forcing.out'
data = np.loadtxt(fname_dpdx)
tarr, ind = np.unique(data[:,0], return_index=True)
if(shuffle_data):
    print('Note: dataset will be randomly shuffled before fixing the time window\n')
    nlen = np.size(data[ind,1][np.where((tarr>tbeg) & (tarr<tend))])
    dpdx_arr = data[ind,1][np.where(tarr>tbeg)]
    np.random.shuffle(dpdx_arr[:])
else:
    dpdx_arr = data[ind,1][np.where((tarr>tbeg) & (tarr<tend))]
    nlen = np.size(dpdx_arr[:])
dpdx = -np.average(dpdx_arr[0:nlen])
#
utau = np.sqrt(dpdx*h/2.)
retau = utau*h/visc
dnu   = visc/utau/h
cf    = utau**2/(ub**2/2.)
np.savetxt(resultsdir+'stats.txt',np.c_[retau,utau,dnu])
#
print("Pressure gradient = ", dpdx)
print("u_tau/u_bulk = ", utau)
print("dnu/h        = ", dnu)
print("Friction Reynolds number = ", retau)
#
# compute single point statistics
#
saves_log = np.loadtxt(datadir+'time.out')
step_log, ind = np.unique(saves_log[:,0], return_index=True)
time_log = saves_log[ind,2]
def find_closest_save(time_log,step_log,time,stp):
    index = np.searchsorted(time_log, time, side="right")
    if(index < np.size(step_log[:])):
        isave = int(step_log[index])-stp//2
        isave = isave - isave%stp
    else:
        isave = -1
    return isave
fldbeg = find_closest_save(time_log,step_log,tbeg,fldstp)
fldend = find_closest_save(time_log,step_log,tend,fldstp)
files = np.sort(glob.glob(datadir+'velstats_fld_???????.out'))
lastfile = files[-1]
fldmax = list(map(int,re.findall(r'([0-9]{7,})', lastfile)))[-1]
if(fldend == -1):
    fldend = fldmax
    index = np.searchsorted(step_log[:], fldmax, side="left")
    print("\nN.B.: Maximum time exceeds number of saves, overwriting end time...")
    tend = time_log[index]
firstfile= files[0]
fldmin = list(map(int,re.findall(r'([0-9]{7,})',firstfile)))[-1]
if(fldmin > fldbeg):
    fldbeg = fldmin
    index = np.searchsorted(step_log[:], fldmin, side="left")
    print("\nN.B.: Minimum larger than the first found save, overwriting  time...")
    tbeg = time_log[index]
if(fldbeg > fldend):
    print("\nError in the selected time interval. Aborting...")
    exit()
t_dns = 1.
tb    = (2.*h)/ub
te    = h/utau
print("Max, min, and interval averaging time in half bulk times (h/ub) = ", tbeg         , tend         , (tend-tbeg)         ,"\n")
print("Max, min, and interval averaging time in turnover times         = ", tbeg*t_dns/te, tend*t_dns/te, (tend-tbeg)*t_dns/te,"\n")
#
fname_beg = datadir+'velstats_fld_'
fname_ext = '.out'

nflds = (fldend-fldbeg)//fldstp + 1
flds = np.array(range(fldbeg,fldend+fldstp,fldstp))
if(shuffle_data):
    np.random.shuffle(flds)
flds = flds[0:nflds]
norm = 1/nflds
grid = np.loadtxt(datadir+'geometry.out')
n1 = int(grid[0,1])
n2 = int(grid[0,2])
#
# single point statistics
#
for ifld in flds[:]:
    #
    fname = fname_beg+str(ifld).zfill(7)+fname_ext
    print(fname)
    data = np.loadtxt(fname)
    if ifld == flds[0]:
        data_avg = data
    else:
        data_avg += data
data_avg = data_avg*norm
#
yc    = np.reshape(data    [:,0 ],(n2,n1),order='C')
zc    = np.reshape(data    [:,1 ],(n2,n1),order='C')
u1    = np.reshape(data_avg[:,2 ],(n2,n1),order='C')
v1    = np.reshape(data_avg[:,3 ],(n2,n1),order='C')
w1    = np.reshape(data_avg[:,4 ],(n2,n1),order='C')
u2    = np.reshape(data_avg[:,5 ],(n2,n1),order='C')
v2    = np.reshape(data_avg[:,6 ],(n2,n1),order='C')
w2    = np.reshape(data_avg[:,7 ],(n2,n1),order='C')
uv    = np.reshape(data_avg[:,8 ],(n2,n1),order='C')
uw    = np.reshape(data_avg[:,9 ],(n2,n1),order='C')
vw    = np.reshape(data_avg[:,10],(n2,n1),order='C')
#
# profile folding function
#
def fold_2d(var, cf='C', isym2=1, isym1=1):
  # support cf='C' and even numbers of cells
  isym = isym2*isym1
  if(abs(isym) != 1): exit()
  n2, n1 = var.shape
  s = 0
  if cf == 'F':
    s = 1
  var[0:n2//2-s:+1, 0:n1//2-s:+1] =      (var[     0:n2//2-s:+1, 0:n1//2-s:+1] + \
                                    isym2*var[n2-1-s:n2//2-1:-1, 0:n1//2-s:+1] + \
                                    isym1*var[     0:n2//2-s:+1, n1-1-s:n1//2-1:-1] + \
                                    isym *var[n2-1-s:n2//2-1:-1, n1-1-s:n1//2-1:-1])*0.25
  var[n2//2-s:n2:+1,  0:n1//2-s:+1] = isym2*var[0:n2//2-s:+1, 0:n1//2-s:+1][::-1, :   ]
  var[ 0:n2//2-s:+1, n1//2-s:n1:+1] = isym1*var[0:n2//2-s:+1, 0:n1//2-s:+1][   :, ::-1]
  var[n2//2-s:n2:+1, n1//2-s:n1:+1] = isym *var[0:n2//2-s:+1, 0:n1//2-s:+1][::-1, ::-1]
  return var
#
u1 = fold_2d(u1,cf='C',isym2=+1, isym1=+1)
v1 = fold_2d(v1,cf='C',isym2=+1, isym1=-1)
w1 = fold_2d(w1,cf='C',isym2=-1, isym1=+1)
u2 = fold_2d(u2,cf='C',isym2=+1, isym1=+1)
v2 = fold_2d(v2,cf='C',isym2=+1, isym1=+1)
w2 = fold_2d(w2,cf='C',isym2=+1, isym1=+1)
uv = fold_2d(uv,cf='C',isym2=+1, isym1=-1)
uw = fold_2d(uw,cf='C',isym2=-1, isym1=+1)
vw = fold_2d(vw,cf='C',isym2=-1, isym1=-1)
#
u2[:] -= u1[:]**2
v2[:] -= v1[:]**2
w2[:] -= w1[:]**2
uv[:] -= u1[:]*v1[:]
uw[:] -= u1[:]*w1[:]
vw[:] -= v1[:]*w1[:]
#
data_avg[:,0 ] = np.reshape(yc,(n1*n2),order='C')
data_avg[:,1 ] = np.reshape(zc,(n1*n2),order='C')
data_avg[:,2 ] = np.reshape(u1,(n1*n2),order='C')
data_avg[:,3 ] = np.reshape(v1,(n1*n2),order='C')
data_avg[:,4 ] = np.reshape(w1,(n1*n2),order='C')
data_avg[:,5 ] = np.reshape(u2,(n1*n2),order='C')
data_avg[:,6 ] = np.reshape(v2,(n1*n2),order='C')
data_avg[:,7 ] = np.reshape(w2,(n1*n2),order='C')
data_avg[:,8 ] = np.reshape(uv,(n1*n2),order='C')
data_avg[:,9 ] = np.reshape(uw,(n1*n2),order='C')
data_avg[:,10] = np.reshape(vw,(n1*n2),order='C')
#
fname = resultsdir + 'stats-single-point-duct-' + casename + fname_ext
with open(fname, 'w') as file:
    file.write(f'ZONE I={n1}, J={n2}, DATAPACKING=POINT\n')
    np.savetxt(file, data_avg, fmt='%16.6e', delimiter='')
    # (e16.7,fortran) is equivalent to (16.6e,python)
#
def interp(n2, n1, y, u, h):
    u1_cl = np.zeros(n2)
    for k in range(n2):
      f = interp1d(y[k,n1//2-3:n1//2], u[k,n1//2-3:n1//2], kind='quadratic', fill_value='extrapolate')
      u1_cl[k] = f(h)
    return u1_cl
#
# bisector statistics
#
n2, n1 = zc.shape
yc_cl = yc[:,n1//2-1]
zc_cl = zc[:,n1//2-1]
u1_cl = interp(n2, n1, yc, u1, h)
v1_cl = interp(n2, n1, yc, v1, h)
w1_cl = interp(n2, n1, yc, w1, h)
u2_cl = interp(n2, n1, yc, u2, h)
v2_cl = interp(n2, n1, yc, v2, h)
w2_cl = interp(n2, n1, yc, w2, h)
uv_cl = interp(n2, n1, yc, uv, h)
uw_cl = interp(n2, n1, yc, uw, h)
vw_cl = interp(n2, n1, yc, vw, h)
fname = resultsdir + 'stats-single-point-duct-centerline-' + casename + fname_ext
with open(fname, 'w') as file:
    np.savetxt(file, np.c_[zc_cl[0:n2//2], \
                           u1_cl[0:n2//2],v1_cl[0:n2//2],w1_cl[0:n2//2], \
                           u2_cl[0:n2//2],v2_cl[0:n2//2],w2_cl[0:n2//2], \
                           uv_cl[0:n2//2],uw_cl[0:n2//2],vw_cl[0:n2//2]], \
               fmt='%16.6e', delimiter='')
#
# diagonal statistics, uniform grid in the cross-section
#
yc_diag = np.diag(yc)
zc_diag = np.diag(zc)
u1_diag = np.diag(u1)
v1_diag = np.diag(v1)
w1_diag = np.diag(w1)
u2_diag = np.diag(u2)
v2_diag = np.diag(v2)
w2_diag = np.diag(w2)
uv_diag = np.diag(uv)
uw_diag = np.diag(uw)
vw_diag = np.diag(vw)
fname = resultsdir + 'stats-single-point-duct-diagonal-' + casename + fname_ext
with open(fname, 'w') as file:
    np.savetxt(file, np.c_[zc_diag[0:n2//2], \
                           u1_diag[0:n2//2],v1_diag[0:n2//2],w1_diag[0:n2//2], \
                           u2_diag[0:n2//2],v2_diag[0:n2//2],w2_diag[0:n2//2], \
                           uv_diag[0:n2//2],uw_diag[0:n2//2],vw_diag[0:n2//2]], \
               fmt='%16.6e', delimiter='')