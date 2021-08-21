"""
Revision date: 2021.08.08

SPDX-License-Identifier: BSD-2-Clause
Copyright (c) 2021 Stuart Nolan. All rights reserved.
"""

import sys
# https://docs.libreoffice.org/pyuno.html
# to debug, import pdb, use pdb.set_trace() and run libreOffice from cmd line
# so that the pdb console will be visible
#import pdb
from scipy.interpolate import interp1d, griddata
from scipy import meshgrid

def psf(method,*args,**kwargs):
    """
    demonstration method for function dispatcher from python module

    will not work from libreOffice as __name__ and __import__ are changed

    To function dispatch from libreOffice, define methods in a class, create
    an instance of that class, and refer to that class instance in getattr

    usage:
        psf("2arg","arg1","arg2") 
        returns: arg1 is 1 and arg2 is argument 2
    """
    func = getattr(__import__(__name__),'m_%s' % method, None)
    res=func(*args,**kwargs)
    return res

def m_2arg(arg1,arg2,arg3=0):
    assert arg1 is not None, "arg1: %s" % arg1
    print('arg1: %s; arg2: %s; arg3: %s' % (arg1,arg2,arg3))

def interp(xVal,xVals,yVals,kind):
    #check xVal is numeric
    #pdb.set_trace()
    kinds = ["linear","cubic","slinear","quadratic","nearest","zero",
             "previous","next"]
    if type(kind) is int or kind in kinds:
        xVals = tuple(el for tup in xVals for el in tup)
        yVals = tuple(el for tup in yVals for el in tup)
        nanMask = [True if type(el) in [int,float] else False for el in yVals]
        yVals = [el for idx, el in enumerate(yVals) if nanMask[idx]]
        xVals = [el for idx, el in enumerate(xVals) if nanMask[idx]]
        if len(xVals) == len(yVals) and len(yVals) > 1:
            f = interp1d(xVals,yVals,kind=kind,fill_value="extrapolate")
            return float(f(xVal))
        else:
            return "insuffcient or invalid data"
    else:
        return "kind: %s is not an integer or in %s" % (kind,kinds)

def interp2d(xvs,yvs,data,xi,yi,meth):
    """
    Wraper function for scipy.interpolate.griddata for use with libreOffice.
    This function will take ranges from libreOffice (passed to python as a
    tuple of tuples), filter out non numeric points in data (including empty 
    cells passed as None), and interpolate the data at one point (xi, yi) via 
    meth.

    Parameters:
        xvs, tuple of "x" data point locations (columns)
        yvs, tuple of "y" data point locations (rows)
        data, data point values at every element of (xvs,yvs)
        xi, x location to interpolate (scaler numeric only)
        yi, y location to interpolate (scaler numeric only)
        meth, griddata interpolation method
    Returns:
        interpolated value from data at point xi, yi

    Ref: https://stackoverflow.com/questions/37662180/interpolate-missing-values-2d-python    
    """
    meths = ["linear","cubic","nearest"]
    if meth in meths:
        #pdb.set_trace()
        xvs = tuple(el for tup in xvs for el in tup)
        yvs = tuple(el for tup in yvs for el in tup)
        data = tuple(el for tup in data for el in tup)
        nanMask = [True if type(el) in [int,float] else False for el in data]
        data = tuple(el for idx,el in enumerate(data) if nanMask[idx])
        xvs,yvs = meshgrid(xvs,yvs)
        xvs = xvs.flatten()[nanMask]
        yvs = yvs.flatten()[nanMask]
        return float(griddata((xvs,yvs),data,(xi,yi),meth))
    else:
        return "meth: %s is not in %s" % (meth,meths)

def pySysVersion():
    """
    """
    return "python path: "+sys.executable+"; Version: "+sys.version.splitlines()[0]

def loArray(tup):
   """
   a range from libre office is passed to python as a 2 diminsiong tuple
   i.e. ((A1,B1,C1),
         (A2,B2,C2)) ~ A1:C2
   
   empty cells in a libreOffice range evaluate to None in python

   empty cells in the transpose of a libreOffice range evaluate to 0.0 in 
   python
   """
   nRows = len(tup)
   nCols = len(tup[0])
   
   return "R: %s; C: %s; last elem: %s;" % (nRows,nCols,tup[nRows-1][nCols-1])
   """
   if isinstance(tup, tuple):
      return "tuple"
   else:
      return "not tuple"
   """
