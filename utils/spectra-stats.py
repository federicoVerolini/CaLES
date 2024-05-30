import numpy as np
precision = 'float64'
use_new_arrays = False
shuffle_data   = False
import sys
args = sys.argv[1:]
casedir = args[3] if len(args) > 3 else ''
datadir = casedir+'data/'
print('Data directory: ' + datadir+'\n')
resultsdir = casedir+'results/'
import os
os.makedirs(resultsdir,exist_ok=True)
#
h = 1.
reb = 39600.
ub = 1.
visc = 2*h*ub/reb
tbeg   = float(args[0]) if len(args) > 0 else 120.
tend   = float(args[1]) if len(args) > 1 else 1.e9
fldstp = int(  args[2]) if len(args) > 2 else 1000
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
import glob
import re
fnamekey = 'spectra_!_fld_???????'
fnamekey_u = fnamekey.replace('!','u',1)
files_spectra_2d_u   = np.sort(glob.glob(datadir+fnamekey_u+'_2d.bin'      ))
files_spectra_meta_u = np.sort(glob.glob(datadir+fnamekey_u+'_kplanes.out' ))
files_spectra_1d_x_u = np.sort(glob.glob(datadir+fnamekey_u+'_1d_x.bin'    ))
files_spectra_1d_y_u = np.sort(glob.glob(datadir+fnamekey_u+'_1d_y.bin'    ))
fnamekey_v = fnamekey.replace('!','v',1)
files_spectra_2d_v   = np.sort(glob.glob(datadir+fnamekey_v+'_2d.bin'      ))
files_spectra_meta_v = np.sort(glob.glob(datadir+fnamekey_v+'_kplanes.out' ))
files_spectra_1d_x_v = np.sort(glob.glob(datadir+fnamekey_v+'_1d_x.bin'    ))
fnamekey_w = fnamekey.replace('!','w',1)
files_spectra_2d_w   = np.sort(glob.glob(datadir+fnamekey_w+'_2d.bin'      ))
files_spectra_meta_w = np.sort(glob.glob(datadir+fnamekey_w+'_kplanes.out' ))
files_spectra_1d_x_w = np.sort(glob.glob(datadir+fnamekey_w+'_1d_x.bin'    ))
fnamekey_p = fnamekey.replace('!','p',1)
files_spectra_2d_p   = np.sort(glob.glob(datadir+fnamekey_p+'_2d.bin'      ))
files_spectra_meta_p = np.sort(glob.glob(datadir+fnamekey_p+'_kplanes.out' ))
files_spectra_1d_x_p = np.sort(glob.glob(datadir+fnamekey_p+'_1d_x.bin'    ))

lastfile = files_spectra_2d_u[-1]
fldmax = list(map(int,re.findall(r'([0-9]{7,})', lastfile)))[-1]
if(fldend == -1):
    fldend = fldmax
    index = np.searchsorted(step_log[:], fldmax, side="left")
    print("\nN.B.: Maximum time exceeds number of saves, overwriting end time...")
    tend = time_log[index]
firstfile= files_spectra_2d_u[0]
fldmin = list(map(int,re.findall(r'([0-9]{7,})',firstfile)))[-1]
if(fldmin > fldbeg):
    fldbeg = fldmin
    index = np.searchsorted(step_log[:], fldmin, side="left")
    print("\nN.B.: Minimum larger than the first found save, overwriting  time...")
    tbeg = time_log[index]
if(fldbeg > fldend):
    print("\nError in the selected time interval. Aborting...")
    exit()
print("Max, min, and interval averaging time in half bulk times (h/ub) =     ", tbeg            , tend            , (tend-tbeg)            ,"\n")
ret = 0.09*reb**0.88
print("Max, min, and interval averaging time in turnover times (estimated) = ", tbeg*ret/(reb/2), tend*ret/(reb/2), (tend-tbeg)*ret/(reb/2),"\n")
casename = '1000'.zfill(5)
#
fname_beg = datadir+'pdfs_fld_'
fname_ext = '.bin'

nflds = (fldend-fldbeg)//fldstp + 1
flds = np.array(range(fldbeg,fldend+fldstp,fldstp))
if(shuffle_data):
    np.random.shuffle(flds)
#
# read grid data
#
flds = flds[0:nflds]
norm = 1/nflds
grid = np.loadtxt(datadir+'geometry.out')
nz   = int(grid[0,2])
f      = open(datadir+'grid.bin','rb')
grid_z = np.fromfile(f,dtype='float64'); f.close()
z_c = np.reshape(grid_z,(nz,4),order='F')[:,2]
z_f = np.reshape(grid_z,(nz,4),order='F')[:,3]
n   = grid[0,:].astype('int')
l   = grid[1,:]
dl  = l[:]/n[:]
#
# read 2d spectra metadata (all files should be the same)
#
data_meta_spectra_u = np.loadtxt(files_spectra_meta_u[0])
nplanes_u = np.size(data_meta_spectra_u[:,0])
#
# spectra
#
def average_spectra_2d(flds,files_2d,n,nplanes,prec='float32'):
    s_2d = np.zeros([n[0]//2+1,n[1]//2+1,nplanes])
    nflds = np.size(flds)
    for ifld in range(nflds):
        fname = re.sub('([0-9]{7,})',str(flds[ifld]).zfill(7),files_2d[ifld])
        data    = np.fromfile(fname,dtype=prec).reshape(np.shape(s_2d),order='F')
        s_2d += data[:,:,:]/nflds
    s_2d = 0.5*(s_2d[:,:,0:nplanes//2:+1] + s_2d[:,:,nplanes-1:nplanes//2-1:-1])/2.
    return s_2d
def average_spectra_1d(flds,files_1d,ioffset_arg,n,nn,prec='float32'):
    s_1d = np.zeros([n//2+1,nn])
    nflds = np.size(flds)
    rsize = 8 if prec == 'float64' else 4
    ioffset = np.size(s_1d)*rsize if ioffset_arg != 0 else 0
    for ifld in range(nflds):
        fname = re.sub('([0-9]{7,})',str(flds[ifld]).zfill(7),files_1d[ifld])
        data    = np.fromfile(fname,dtype=prec,count=np.size(s_1d),offset=ioffset).reshape(np.shape(s_1d),order='F')
        s_1d += data[:,:]/nflds
    s_1d = 0.5*(s_1d[:,0:nn//2:+1] + s_1d[:,nn-1:nn//2-1:-1])/2.
    return s_1d
def compute_scales(s_1d,dl):
    nn      = np.size(s_1d[:,:],axis=0)
    rr      = np.fft.irfft(s_1d,axis=0)[0:nn]/2.
    var     = rr[0,:]
    rr     /= var[:]
    lint    = np.sum(rr*dl,axis=0)
    d2rdl2  = (rr[1,:]-2*rr[0,:]+rr[1,:])/dl**2
    llambda = (-0.5*d2rdl2[:])**(-.5)
    return rr,lint,llambda
#
s_2d_u    = average_spectra_2d(flds,files_spectra_2d_u  ,  n   ,nplanes_u,prec=precision)
s_1d_x0_u = average_spectra_1d(flds,files_spectra_1d_x_u,0,n[0],nz       ,prec=precision)
s_1d_xs_u = average_spectra_1d(flds,files_spectra_1d_x_u,1,n[0],nz       ,prec=precision)
s_1d_y0_u = average_spectra_1d(flds,files_spectra_1d_y_u,0,n[1],nz       ,prec=precision)
s_1d_ys_u = average_spectra_1d(flds,files_spectra_1d_y_u,1,n[1],nz       ,prec=precision)
r_x_u,lint_x_u,llambda_x_u = compute_scales(s_1d_x0_u,dl[0])
r_y_u,lint_y_u,llambda_y_u = compute_scales(s_1d_y0_u,dl[1])
#
# save data
#
# (1) some grid info
#
kx = np.arange(0,n[0]/2.+1,1)*2*np.pi/l[0]       # wavenumbers
ky = np.arange(0,n[1]/2.+1,1)*2*np.pi/l[1]
with np.errstate(divide='ignore'):
    lambx = 2.*np.pi/kx[:];lambx[0] = lambx[1]   # wavelengths
    lamby = 2.*np.pi/ky[:];lamby[0] = lamby[1]
rx = np.arange(0.,l[0]/2.+dl[0],dl[0])           # separations
ry = np.arange(0.,l[1]/2.+dl[1],dl[1])
fname = resultsdir + 'stats-spectra-chan-' + casename + '-xgrids' + '.out'
np.savetxt(fname,np.c_[kx,lambx,rx])
fname = resultsdir + 'stats-spectra-chan-' + casename + '-ygrids' + '.out'
np.savetxt(fname,np.c_[ky,lamby,ry])
fname = resultsdir + 'stats-spectra-chan-' + casename + '-zgrids' + '.out'
np.savetxt(fname,np.c_[z_c,z_f])
#
# (2) 1D spectra
#
fname = resultsdir + 'stats-spectra-chan-' + casename + '1d-x0-u' + '.out'
np.savetxt(fname,s_1d_x0_u)
fname = resultsdir + 'stats-spectra-chan-' + casename + '1d-xs-u' + '.out'
np.savetxt(fname,s_1d_xs_u)
fname = resultsdir + 'stats-spectra-chan-' + casename + '1d-y0-u' + '.out'
np.savetxt(fname,s_1d_y0_u)
fname = resultsdir + 'stats-spectra-chan-' + casename + '1d-ys-u' + '.out'
np.savetxt(fname,s_1d_ys_u)
#
# (3) autocorrelations
#
fname = resultsdir + 'stats-spectra-chan-' + casename + '1d-ruu-x0-u' + '.out'
np.savetxt(fname,r_x_u)
fname = resultsdir + 'stats-spectra-chan-' + casename + '1d-ruu-y0-u' + '.out'
np.savetxt(fname,r_y_u)
#
# (4) integral scale and Taylor microscale profiles
#
fname = resultsdir + 'stats-spectra-chan-' + casename + '-scales-x0-u' + '.out'
np.savetxt(fname,np.c_[z_c[0:nz//2],lint_x_u,llambda_x_u])
fname = resultsdir + 'stats-spectra-chan-' + casename + '-scales-y0-u' + '.out'
np.savetxt(fname,np.c_[z_c[0:nz//2],lint_y_u,llambda_y_u])
#
# (5) 2D spectra along planes
#
for k in range(nplanes_u//2):
    kplane = int(data_meta_spectra_u[k,0])
    fname = resultsdir + 'stats-spectra-chan-' + casename + '-2d-u-plane-' + str(kplane).zfill(5) + '.out'
    np.savetxt(fname,s_2d_u[:,:,k])
##
## sandbox below, to be removed
##
#import matplotlib
#import matplotlib.pyplot as plt
#for i in range(n[0]//2+1):
#    for j in range(n[1]//2+1):
#        s_2d_u[i,j,:] = s_2d_u[i,j,:]*kx[i]*ky[j]
#for k in range(nz//2):
#    s_1d_y0_u[:,k] = s_1d_y0_u[:,k]*ky[:]
#plt.xlim(ky[1]/10.,ky[-1])
#plt.ylim(kx[1]/10.,kx[-1])
#plt.xscale('log');plt.yscale('log');plt.contourf(ky,kx,s_2d_u[:,:,2]);plt.title('2D');plt.colorbar();plt.show()
#plt.contourf(kx,z_c[2:nz//2]*ret,np.transpose(s_1d_x0_u[:,2:]));plt.title('1D, x0');plt.colorbar();plt.show()
#plt.contourf(kx,z_c[2:nz//2]*ret,np.transpose(s_1d_xs_u[:,2:]));plt.title('1D, xs');plt.colorbar();plt.show()
##plt.xlim(ry[1]/2.,ry[-1])
#plt.contourf(ry*180.,z_c[2:nz//2]*ret,np.transpose(s_1d_y0_u[:,2:nz//2]));plt.title('1D, y0');plt.colorbar();plt.show()
#plt.contourf(ky,z_c[2:nz//2]*ret,np.transpose(s_1d_ys_u[:,2:]));plt.title('1D, ys');plt.colorbar();plt.show()
##
#plt.contourf(rx*ret,z_c[2:nz//2]*ret,np.transpose(r_x_u[:,2:]),levels=50);plt.colorbar();plt.show()
#plt.contourf(ry*ret,z_c[2:nz//2]*ret,np.transpose(r_y_u[:,2:]),levels=50);plt.colorbar();plt.show()
#plt.plot(lint_x_u*ret,z_c[0:nz//2]*ret);plt.show()
#plt.plot(llambda_x_u*ret,z_c[0:nz//2]);plt.show()
#plt.plot(lint_y_u*ret,z_c[0:nz//2]*ret);plt.show()
#plt.plot(llambda_y_u*ret,z_c[0:nz//2]);plt.show()
