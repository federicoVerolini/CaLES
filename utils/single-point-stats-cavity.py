# -
#
# SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
# SPDX-License-Identifier: MIT
#
# -
#!/usr/bin/env python
def read_single_field_binary(filenamei,iskip):
    import numpy as np
    #
    # setting up some parameters
    #
    iprecision = 8            # precision of the real-valued data
    r0 = np.array([0.,0.,0.]) # domain origin
    non_uniform_grid = True
    precision  = 'float64'
    if(iprecision == 4): precision = 'float32'
    #
    # read geometry file
    #
    geofile  = "geometry.out"
    geo = np.loadtxt(geofile, comments = "!", max_rows = 2)
    ng = geo[0,:].astype('int')
    l  = geo[1,:]
    dl = l/(1.*ng)
    #
    # read and generate grid
    #
    xp = np.arange(r0[0]+dl[0]/2.,r0[0]+l[0],dl[0]) # centered  x grid
    yp = np.arange(r0[1]+dl[1]/2.,r0[1]+l[1],dl[1]) # centered  y grid
    zp = np.arange(r0[2]+dl[2]/2.,r0[2]+l[2],dl[2]) # centered  z grid
    xu = xp + dl[0]/2.                              # staggered x grid
    yv = yp + dl[1]/2.                              # staggered y grid
    zw = zp + dl[2]/2.                              # staggered z grid
    if(non_uniform_grid):
        f   = open('grid.bin','rb')
        grid_z = np.fromfile(f,dtype=precision)
        f.close()
        grid_z = np.reshape(grid_z,(ng[2],4),order='F')
        zp = r0[2] + np.transpose(grid_z[:,2]) # centered  z grid
        zw = r0[2] + np.transpose(grid_z[:,3]) # staggered z grid
    #
    # read binary file
    #
    n           = (ng[:]/iskip[:]).astype(int)
    data        = np.zeros([n[0],n[1],n[2]])
    fld         = np.fromfile(filenamei,dtype=precision)
    data[:,:,:] = np.reshape(fld,(n[0],n[1],n[2]),order='F')
    #
    # reshape grid
    #
    xp = xp[0:ng[0]:iskip[0]]
    yp = yp[0:ng[1]:iskip[1]]
    zp = zp[0:ng[2]:iskip[2]]
    xu = xu[0:ng[0]:iskip[0]]
    yv = yv[0:ng[1]:iskip[1]]
    zw = zw[0:ng[2]:iskip[2]]
    return data, xp, yp, zp, xu, yv, zw

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import numpy as np
    import sys
    import os

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
    #
    fname_ext = '.out'
    iskipx      = 1
    iskipy      = 1
    iskipz      = 1
    iskip       = np.array([iskipx,iskipy,iskipz]).astype(int)
    #
    # u along z direction
    #
    data, xp, yp, zp, xu, yv, zw = read_single_field_binary('vex_fld_0010000.bin',iskip)
    fname = resultsdir + 'stats-single-point-cavi-vertical-' + casename + fname_ext
    np.savetxt(fname, np.column_stack((zp, 0.5*(data[63,63,:]+data[63,64,:])/2.0)), fmt='%16.6e', delimiter='')
    #
    # w along x direction
    #
    data, xp, yp, zp, xu, yv, zw = read_single_field_binary('vez_fld_0010000.bin',iskip)
    fname = resultsdir + 'stats-single-point-cavi-horizontal-' + casename + fname_ext
    np.savetxt(fname, np.column_stack((xp, 0.5*(data[:,63,63]+data[:,64,63])/2.0)), fmt='%16.6e', delimiter='')