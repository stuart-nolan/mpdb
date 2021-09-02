"""
misc_cython.pyx: equations for thermophycial property estimation

To build:

python misc_cython-setup.py build_ext --inplace

Revision Date: 2021.09.02

SPDX-License-Identifier: BSD-2-Clause
Copyright (c) 2021 Stuart Nolan. All rights reserved.
"""

from libc.math cimport log10, exp

cimport cython
@cython.boundscheck(False) # turn off bounds-checking
@cython.wraparound(False)  # turn off negative index wrapping

cpdef enthalpy_Shomate_bounded(float T, list TBounds, list sPs,
                               float TRef=298.15):
    """
    enthalpy_Shomate_bounded(float T, list Tbounds, list sPs, 
                             float TRef=298.15)

    Shomate Equation for enthalpy using multiple temperature bounded 
    parameter sets.  The first parameter set is used for witch 
    Tmin <= T <= Tmax where Tmin and Tmax are the minimum and maximum 
    temperatures over which the parameter set is valid. 

    Parameters
      T, temperature in Kelvin
      TBounds, [[Tmin_0, Tmax_0], [Tmin_1, Tmax_1], [Tmin_2, ...] ... ]
               An ordered list of 2 item lists of minimum and maximum 
               temperatures.  The number of lists in TBounds must equal
               the number of lists in sPs (described below).  List item 
               0 in TBounds corresponds to list item 0 in sPs.  List item 
               1 in TBounds corresponds to list item 1 in sPs, etc... 
      sPs, [[A_0, B_0, C_0, D_0, E_0, F_0, G_0, H_0], [A_1, B_1, C_1, D_1, 
            E_1, F_1, G_1, H_1], [A_2, ... ], ... ]
           An ordered list of 8 item lists of parameters for the
           shomate enthalpy equation (see shomate_enthalpy).  The
           number of lists in sPs must equal the number of lists in
           TBounds (described above).  List item 0 in sPs corresponds
           to list item 0 in TBounds, List item 1 in sPs corresponds
           to list item 1 in TBounds, etc...
      TRef, reference temperature in Kelvin

   Returns
      enthalpy_Shomate(T, sP[index]) - HRef if T and TRef are within TBounds 
      (TRef = 298.15 K is always ok) otherwise None
    """
    cdef float HRef = 0
    cdef int TRef_bounded = 1
    if TRef != 298.15:
        TRef_bounded = 0
        for idx, TB in enumerate(TBounds):
            if TRef >= TB[0] and TRef <= TB[1]:
                HRef = enthalpy_Shomate(TRef,sPs[idx])
                TRef_bounded = 1

        if not TRef_bounded:
            return None

    for idx, TB in enumerate(TBounds):
        if T >= TB[0] and T <= TB[1]:
            return enthalpy_Shomate(T, sPs[idx]) - HRef

    return None

cpdef float enthalpy_Shomate(float T, list sP):
    """
    enthalpy_Shomate(float T, list sP)
    
    Enthalpy at T relative to enthalpy at 298.15 from Nist Shomate equation
      t=T/1000
      enthalpy_Shomate = A*t + B/2*t^2 + C/3*t^3 + D/4*T^4 - E*1/t + F - H
    Parameters
      T, temperature in Kelvin
      sP, A=sP[0], B=sP[1], C=sP[2], D=sP[3], E=sP[4]
      A, B, C, D, E, F, & H are regression coefficients
    Returns
      enthalpy change in kJ/mol between T and TRef = 298.15 K
    """
    cdef float t = T/1000
    return sP[0]*t + sP[1]/2.0*t**2 + sP[2]/3.0*t**3 + sP[3]/4.0*t**4 - sP[4]/t + sP[5] - sP[7]

cpdef float enthalpy_f1(float T, list hCP, float TRef=298.15):
    """
    enthalpy_f1(float T, list hCP, float TRef=298.15)

    enthalpy change correlation (form 1: slop vapor & solid data)
      h_f1 (J/mol) = A*(T-TRef) + 1/2*B*(T^2-TRef^2) - C*(1/T-1/TRef)
    
    Parameters
      T, temperature in Kelvin
      TRef, reference temperature in Kelvin
      hCP, A=hCP[0], B=hCP[1], C=hCP[2]
      A, B, and C are regression coefficients  

    Returns
      enthalpy change in J/mol/K between T and TRef
    """
    return hCP[0]*(T-TRef) + 1.0/2.0*hCP[1]*(T**2-TRef**2) - hCP[2]*(1/T-1/TRef)

cpdef float heatCapacity_f1(float T, list hCP):
    """
    heatCapacity_f1(float T, list hCP) (slop vapors & solids phases)

    heatCapacity_f1 (J/mol/K) = A + B*T + C/T^2
    
    Parameters
      T, temperature in Kelvin
      hCP, A=hCP[0], B=hCP[1], C=hCP[2]
      A, B, and C are regression coefficients  

    Returns
      heat capacity in J/mol/K at T
    """
    return hCP[0] + hCP[1]*T + hCP[2]/T**2

cpdef float heatCapacity_Shomate(float T, list hCP):
    """
    heatCapacity_Shomate(float T, list hCP)
      NIST vapor, liquid, and solid phases

    heat capacity correlation
      hC (J/mol/K) = A + B*t + C*t^2 + D*t^3 + E/t^2
      t (K) = T/1000.0
    
    Parameters
      T, temperature in Kelvin
      hCP, A=hCP[0], B=hCP[1], C=hCP[2], D=hCP[3], E=hCP[4]
      A, B, C, D, and E are regression coefficients  

    Returns
      heat capacity in J/mol/K at T
    """
    cdef float t = T/1000.0 
    return hCP[0] + hCP[1]*t + hCP[2]*t**2 + hCP[3]*t**3 + hCP[4]/t**2

cpdef float henry_rSSI(float T, list hbpSIP):
    """
    henry_rSSI(float T, list hbpSIP)

    Henry's law correlation in molality per pressure SI units from R. Sander
    
    henry_rSSI (mol/(kg solvent)/Pa) = Ho*exp({dln(H)/d(1/T)}*(1/T-1/To))
      
      where
        dln(H)/d(1/T) = -deltaH_sol/R,
        To = reference tempeature (298.15 K)
    
    Parameters
      T, temperature in Kelvin
      hbpSIP, Ho=hbpSIP[0], dln(H)/d(1/T)=hbpSIP[1]
      Ho and dln(H)/d(1/T) are regression coefficients

    Returns
      Henry's law constant at T in mol/kg/Pa
    """
    return hbpSIP[0]*exp(hbpSIP[1]*(1.0/T-1.0/298.15))

cpdef float mCVTC(list molFrac, list vTC, list vV, list fW, float eps=1.0):
    """
    mCVTC(list molFrac, list vTC, list vV, list fW, float eps=1.0)
    
    multi-Component Vapor Thermal Conductivity in W/m/K

    Parameters:
        molFrac, list of mol fractions
        vTC, list of pure component vapor thermal conductivities (components 
             ordered as in molFrac, fW, and vV) 
        vV, list of pure component vapor viscosities (components ordered as
            in molFrac, fW, and vTC)
        fW, list of formula weights for each component in molFrac
        eps, epsilon - typical value is 1.0
    
    Returns:
        multi-component vapor thermal conductivity in W/m/K
    """
    cdef list Aij = []
    cdef float tCm = 0, fWR, sumyjAij, vR
    cdef int idx, jdx, nMF
    
    nMF = len(molFrac)
    
    for idx in range(nMF):
        Aij.append([])
        sumyjAij = 0
        for jdx in range(nMF):
            fWR=fW[idx]/fW[jdx]
            vR = vV[idx]/vV[jdx]
            if idx == jdx:
                sumyjAij += molFrac[jdx]
            elif idx < jdx:
                Aij[idx].append(eps*(1+vR**(1./2)*(1/fWR)**(1./4))**2/(8*(1+fWR))**(1./2))
                sumyjAij += molFrac[jdx]*Aij[idx][jdx-(1+idx)]
            else:
                sumyjAij += molFrac[jdx]*Aij[jdx][idx-(1+jdx)]*vR/fWR
        tCm += molFrac[idx]*vTC[idx]/sumyjAij
    
    return tCm

cpdef float mCVV(list molFrac, list vV, list fW):
    """
    mCVV(list molFrac, list vV, list fW)
    
    multi-Component Vapor Viscosity in micropoise
    
    Parameters:
        molFrac, list of mol fractions
        vV, list of pure component vapor viscosities (components ordered as
            in fW mol)
        fW, list of formula weights for each component in molFrac
    
    Returns:
        multi-component vapor viscosity in micropoise
    """
    cdef list Aij = []
    cdef float vm = 0, fWR, sumyjAij, vR
    cdef int idx, jdx, nMF
    
    nMF = len(molFrac)
    
    for idx in range(nMF): #row
        Aij.append([])
        sumyjAij = 0
        for jdx in range(nMF): #col
            fWR=fW[idx]/fW[jdx]
            vR = vV[idx]/vV[jdx]
            if idx == jdx:
                sumyjAij += molFrac[jdx]
            elif idx < jdx:
                Aij[idx].append((1+vR**(1./2)*(1/fWR)**(1./4))**2/(8*(1+fWR))**(1./2))
                sumyjAij += molFrac[jdx]*Aij[idx][jdx-(1+idx)]
            else:
                sumyjAij += molFrac[jdx]*Aij[jdx][idx-(1+jdx)]*vR/fWR
        vm += molFrac[idx]*vV[idx]/sumyjAij
    
    return vm

cpdef float mFWQ(list mFL, list qL):
    """
    mFWQ(list mFL, list qL):
        
    (mole or mass) Fraction Weighted Quantity
    
    Parameters:
        mFL, list of mole or mass fractions, sum(mFL) = 1
        qL, list of quantities corresponding to items in mFL
    
    Returns:
        weighted averaged of items in qL
    """
    cdef float mF, aveQ = 0
    cdef int idx
    
    for idx,mF in enumerate(mFL):
        aveQ += mF*qL[idx]
    
    return aveQ
