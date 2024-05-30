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
fldstp = int(  args[2]) if len(args) > 2 else 100
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
fnamekey = 'pdfs_fld_???????'
files_meta_pdf  = np.sort(glob.glob(datadir+fnamekey+'_pdf_meta.out' ))
files_pdf       = np.sort(glob.glob(datadir+fnamekey+'_pdf.bin'      ))
files_meta_jpdf = np.sort(glob.glob(datadir+fnamekey+'_jpdf_meta.out'))
files_jpdf      = np.sort(glob.glob(datadir+fnamekey+'_jpdf.bin'     ))
lastfile = files_pdf[-1]
fldmax = list(map(int,re.findall(r'([0-9]{7,})', lastfile)))[-1]
if(fldend == -1):
    fldend = fldmax
    index = np.searchsorted(step_log[:], fldmax, side="left")
    print("\nN.B.: Maximum time exceeds number of saves, overwriting end time...")
    tend = time_log[index]
firstfile= files_pdf[0]
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
z_g = np.reshape(grid_z,(nz,4),order='F')[:,2]
#
# read pdf metadata (all files should be the same)
#
data_meta_pdf = np.loadtxt(files_meta_pdf[0])
npdf   = int(data_meta_pdf[0])
nvars  = int(np.size(data_meta_pdf[:]-1)/2)
pdfmin = data_meta_pdf[1:(nvars+1)*2:2]
pdfmax = data_meta_pdf[2:(nvars+1)*2:2]
dpdf   = (pdfmax-pdfmin)/npdf
#
# pdfs
#
iu   = 0
iv   = 1
iw   = 2
ix   = 3
ip   = 4
#
u   = np.zeros([npdf,nz])
v   = np.zeros([npdf,nz])
w   = np.zeros([npdf,nz])
p   = np.zeros([npdf,nz])
#
for ifld in flds[:]:
    #
    fname = fname_beg+str(ifld).zfill(7)+'_pdf'+fname_ext
    rawdata = np.fromfile(fname,dtype=precision)
    data    = np.reshape(rawdata,(npdf,nz,nvars),order='F')
    u[:,:] += data[:,:,iu]
    v[:,:] += data[:,:,iv]
    w[:,:] += data[:,:,iw]
    p[:,:] += data[:,:,ip]
#
# merging mirror-symmetry planes
#
u = u[:,0:nz//2:+1] + u[:,nz-1:nz//2-1:-1]
v = v[:,0:nz//2:+1] + v[:,nz-1:nz//2-1:-1]
w = w[:,0:nz//2:+1] + w[:,nz-1:nz//2-1:-1]
p = p[:,0:nz//2:+1] + p[:,nz-1:nz//2-1:-1]
def pdf(var,dpdf):
    norm = np.sum(var[:,:],axis=0)
    norm[np.where(norm == 0.)] = 1.
    var[:,:] /= norm[:]*dpdf
    return var
u[:,:] = pdf(u,dpdf[iu])
v[:,:] = pdf(v,dpdf[iv])
w[:,:] = pdf(w,dpdf[iw])
p[:,:] = pdf(p,dpdf[ip])
#
xu = np.arange(pdfmin[iu]+dpdf[iu]/2,pdfmax[iu],dpdf[iu]) # to be done in a separate script, for plotting
xv = np.arange(pdfmin[iv]+dpdf[iv]/2,pdfmax[iv],dpdf[iv]) # to be done in a separate script, for plotting
xw = np.arange(pdfmin[iw]+dpdf[iw]/2,pdfmax[iw],dpdf[iw]) # to be done in a separate script, for plotting
xp = np.arange(pdfmin[ip]+dpdf[ip]/2,pdfmax[ip],dpdf[ip]) # to be done in a separate script, for plotting
fname = resultsdir + 'stats-pdf-chan-u-' + casename + '.out'
np.savetxt(fname,u)
fname = resultsdir + 'stats-pdf-chan-v-' + casename + '.out'
np.savetxt(fname,v)
fname = resultsdir + 'stats-pdf-chan-w-' + casename + '.out'
np.savetxt(fname,w)
fname = resultsdir + 'stats-pdf-chan-p-' + casename + '.out'
np.savetxt(fname,p)
#
# read joint pdf metadata (all files should be the same)
#
data_meta_jpdf = np.loadtxt(files_meta_jpdf[0])
njpdf   = int(data_meta_jpdf[:,0][0])
nvars   = int((np.shape(data_meta_jpdf[:,:])[1]-3)/2)
nplanes = int((np.shape(data_meta_jpdf[:,:])[0]))
pdfmin = data_meta_jpdf[:,3:(nvars+3)*2:2]
pdfmax = data_meta_jpdf[:,4:(nvars+3)*2:2]
dpdf   = (pdfmax-pdfmin)/npdf
#
# pdfs
#
uv  = np.zeros([njpdf,njpdf,nplanes])
uw  = np.zeros([njpdf,njpdf,nplanes])
up  = np.zeros([njpdf,njpdf,nplanes])
vw  = np.zeros([njpdf,njpdf,nplanes])
vp  = np.zeros([njpdf,njpdf,nplanes])
wp  = np.zeros([njpdf,njpdf,nplanes])
#
for ifld in flds[:]:
    #
    fname = fname_beg+str(ifld).zfill(7)+'_jpdf'+fname_ext
    rawdata = np.fromfile(fname,dtype=precision)
    data    = np.reshape(rawdata,(npdf,npdf,nplanes,nvars,nvars),order='F')
    uv[:,:,:] += data[:,:,:,iu,iv]
    uw[:,:,:] += data[:,:,:,iu,iw]
    up[:,:,:] += data[:,:,:,iu,ip]
    vw[:,:,:] += data[:,:,:,iv,iw]
    vp[:,:,:] += data[:,:,:,iv,ip]
    wp[:,:,:] += data[:,:,:,iw,ip]
def pdf(var,djpdf):
    norm = np.sum(var[:,:,:],axis=(0,1))
    norm[np.where(norm == 0.)] = 1.
    var[:,:,:] /= norm[:]*djpdf[:]
    return var
djpdf_u = dpdf[:,iu]
djpdf_v = dpdf[:,iv]
djpdf_w = dpdf[:,iw]
djpdf_p = dpdf[:,ip]
uv[:,:,:] = pdf(uv,djpdf_u*djpdf_v)
uw[:,:,:] = pdf(uw,djpdf_u*djpdf_w)
up[:,:,:] = pdf(up,djpdf_u*djpdf_p)
vw[:,:,:] = pdf(vw,djpdf_v*djpdf_w)
vp[:,:,:] = pdf(vp,djpdf_v*djpdf_p)
wp[:,:,:] = pdf(wp,djpdf_w*djpdf_p)
#
# n.b.: when plotting xu, xv, xw, xp, we can subtract the local mean values, fetched from the single-point stats
#
xu = np.arange(pdfmin[0,iu]+dpdf[0,iu]/2,pdfmax[0,iu],dpdf[0,iu]) # to be done in a separate script, for plotting
xv = np.arange(pdfmin[0,iv]+dpdf[0,iv]/2,pdfmax[0,iv],dpdf[0,iv]) # to be done in a separate script, for plotting
xw = np.arange(pdfmin[0,iw]+dpdf[0,iw]/2,pdfmax[0,iw],dpdf[0,iw]) # to be done in a separate script, for plotting
xp = np.arange(pdfmin[0,ip]+dpdf[0,ip]/2,pdfmax[0,ip],dpdf[0,ip]) # to be done in a separate script, for plotting
#
# save data plane-by-plane
#
for k in range(nplanes):
    kplane = int(data_meta_jpdf[k,1])
    fname = resultsdir + 'stats-jpdf-chan-' + casename + '-uv-plane-' + str(kplane).zfill(5) + '.out'
    np.savetxt(fname,uv[:,:,k])
    fname = resultsdir + 'stats-jpdf-chan-' + casename + '-uw-plane-' + str(kplane).zfill(5) + '.out'
    np.savetxt(fname,uw[:,:,k])
    fname = resultsdir + 'stats-jpdf-chan-' + casename + '-up-plane-' + str(kplane).zfill(5) + '.out'
    np.savetxt(fname,up[:,:,k])
    fname = resultsdir + 'stats-jpdf-chan-' + casename + '-vw-plane-' + str(kplane).zfill(5) + '.out'
    np.savetxt(fname,vw[:,:,k])
    fname = resultsdir + 'stats-jpdf-chan-' + casename + '-vp-plane-' + str(kplane).zfill(5) + '.out'
    np.savetxt(fname,vp[:,:,k])
    fname = resultsdir + 'stats-jpdf-chan-' + casename + '-wp-plane-' + str(kplane).zfill(5) + '.out'
    np.savetxt(fname,wp[:,:,k])
