
#------------------------------------------------------------------------------------------------------------
# function [ table ] = subs_getresults(ire, icp, iblk, icol)
#   or
# function [ blk ] = subs_getresults(ire, icp, iblk, 'all')
#
# return selected columns of the results of the subsurface stress calculation for one block for a
# contact problem.
#
#  iblk              - number of elements in potential contact area
#  icol(ncol)        - requested columns of result data
#                       1-- 3: x,y,z        - positions of points for subsurface stress calculation
#                       4-- 6: ux,uy,uz     - elastic displacements
#                       7-- 9: sighyd,vm,tr - hydrostatic, von Mises and Tresca stresses
#                      10--12: sigma1,2,3   - principal stresses
#                      13--15: sigxx,yx,zx  - components of the full stress tensor
#                      16--18: sigxy,yy,zy
#                      19--21: sigxz,yz,zz
#  table(npnt,ncol)  - result data, with npnt = nx*ny*nz entries filled per column
#  blk               - structure with subsurface results, as used by loadstrs / plotstrs
#------------------------------------------------------------------------------------------------------------

# Copyright 2008-2023 by Vtech CMCC.
#
# Licensed under Apache License v2.0.  See the file "LICENSE.txt" for more information.

# category 6: m=*, cp     - require icp>0, default 1

import sys
import numpy               as     np
from python_intfc          import cntc_dll
from ctypes                import c_int, c_double, POINTER
from .subs_getblocksize    import subs_getblocksize

def subs_getresults(ire, icp, iblk, icol=None):

    totcol = 21

    if (not isinstance(ire, int)):
        ire = 1
    if (not isinstance(icp, int)):
        icp = 1
    if (icp<=0):
        sys.exit('ERROR in subs_getresults: not available for icp=%d' % icp)
    if (not isinstance(iblk, int)):
        iblk = 1

    # convert icol to NumPy ndarray, allowing for scalar or list of scalars
    if (not isinstance(icol, np.ndarray)):
        icol = np.atleast_1d( np.array( icol, dtype=c_int ) )

    # TODO: support 'all' option

  # if (nargin>=4 & ischar(icol))
  #    make_struct = 1;
  #    icol = [1:totcol];

    if (np.any(icol<=0) or np.any(icol>totcol)):
        print('icol=', icol)
        sys.exit('ERROR in subs_getresults: columns must be 1 <= icol <= %d' % totcol)

    # get block-size
    nx, ny, nz = subs_getblocksize(ire, icp, iblk)

    if (nx<=0 or ny<=0 or nz<=0):
        print('nx,ny,nz=', [nx, ny, nz])
        sys.exit('ERROR in subs_getresults: no data for iblk = %d' % iblk)


    make_struct = 0

    npnt = nx * ny * nz
    ncol = len(icol)

    # TODO: check row/column storage, transpose/reshape?

    table = np.zeros( (ncol,npnt), dtype=c_double)

    cntc_dll.subs_getresults(c_int(ire), c_int(icp), c_int(iblk), c_int(npnt), c_int(ncol),
                                                                  icol.ctypes.data_as(POINTER(c_int)),
                                                                  table.ctypes.data_as(POINTER(c_double)))

    # TODO: conversion to struct cf. plotstrs

  # if (make_struct)
  #     blk = struct('nx',nx, 'ny',ny, 'nz',nz);
  
  #     # sort data such that x runs fastest, then y, then z
  #     [~, iperm] = sort(table(:,1)); table = table(iperm,:);
  #     [~, iperm] = sort(table(:,2)); table = table(iperm,:);
  #     [~, iperm] = sort(table(:,3)); table = table(iperm,:);
  
  #     # the data-lines are given with iz running fastest, ix slowest.
  #     blk.npoints = nx * ny * nz;
  #     blk.x       = reshape(table(:, 1), nx, ny, nz); blk.x = squeeze(blk.x(:,1,1));
  #     blk.y       = reshape(table(:, 2), nx, ny, nz); blk.y = squeeze(blk.y(1,:,1)); blk.y = blk.y';
  #     blk.z       = reshape(table(:, 3), nx, ny, nz); blk.z = squeeze(blk.z(1,1,:));
  #     blk.ux      = reshape(table(:, 4), nx, ny, nz);
  #     blk.uy      = reshape(table(:, 5), nx, ny, nz);
  #     blk.uz      = reshape(table(:, 6), nx, ny, nz);
  #     blk.sighyd  = reshape(table(:, 7), nx, ny, nz);
  #     blk.sigvm   = reshape(table(:, 8), nx, ny, nz);
  #     blk.sigtr   = reshape(table(:, 9), nx, ny, nz);
  #     blk.sigma1  = reshape(table(:,10), nx, ny, nz);
  #     blk.sigma2  = reshape(table(:,11), nx, ny, nz);
  #     blk.sigma3  = reshape(table(:,12), nx, ny, nz);
  #     blk.sigxx   = reshape(table(:,13), nx, ny, nz);
  #     blk.sigxy   = reshape(table(:,16), nx, ny, nz);
  #     blk.sigyy   = reshape(table(:,17), nx, ny, nz);
  #     blk.sigxz   = reshape(table(:,19), nx, ny, nz);
  #     blk.sigyz   = reshape(table(:,20), nx, ny, nz);
  #     blk.sigzz   = reshape(table(:,21), nx, ny, nz);
  #     table = blk;
  
  # end # make_struct

    return table

# end function subs_getresults

#------------------------------------------------------------------------------------------------------------

