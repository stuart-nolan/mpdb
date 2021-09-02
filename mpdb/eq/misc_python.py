"""
mpdb: Material Property Data Base as python module
misc_python.py: utility methods for thermophysical property estimation

Revision Date: 2021.09.02

SPDX-License-Identifier: BSD-2-Clause
Copyright (c) 2021 Stuart Nolan. All rights reserved.
"""
import CoolProp as cP
from mpdb.eq import yaws_python as yp

def antoine(T, vPP):
    """
    antoine(T, vPP)
      vapor pressure units are dependent on units of T and vPP
      NIST, T in K, vP in bar
      DECHEMA, T in deg. C, vP in mmHg

      Antoine vapor pressure = 10**(A - (B/(T + C)))
    
    Parameters
      T, temperature
      vPP, A=vPP[0], B=vPP[1], C=vPP[2]
      A, B, and C are regression coefficients  

    Returns
      vapor pressure at T
    """
    return 10**(vPP[0] - (vPP[1]/(T+vPP[2])))

def coolPropFluids():
    """
    Returns a list of tuples [('coolprop fluid CAS', 'coolprop fluid name'),
                              ...]
    for all the fluids in available in CoolProp
    alt function to get coolprop fluids:
        cPFluids = cP.CoolProp.FluidsList()
    """
    return dict([(cP.CoolProp.get_fluid_param_string(fluid,"CAS"), fluid)
                 for fluid in 
                 cP.CoolProp.get_global_param_string("fluids_list").split(',')])
    
class coolPropEqns:
    #cP = __import__('CoolProp')
    def __init__(self,component):
        """
        Thermophysical equations similar to yaws equations in 
        mpdb.eq.yaws_python as calculated by CoolProp

        NOTE: CoolProp does not implement all thermophysical quantites.
              E.g. IsoButene viscosity:
              import CoolProp.CoolProp as cP
              cP.PropsSI('V','T',298.15,'P',101325,'IsoButene')
        
        from CoolProp/include/DataStructures.h
        0 iphase_liquid, < Subcritical liquid
        1 iphase_supercritical, < Supercritical (p > pc, T > Tc)
        2 iphase_supercritical_gas, < Supercritical gas (p < pc, T > Tc)
        3 iphase_supercritical_liquid, < Supercritical liquid (p > pc, T < Tc)
        4 iphase_critical_point, < At the critical point
        5 iphase_gas, < Subcritical gas
        6 iphase_twophase, < Twophase
        7 iphase_unknown, < Unknown phase
        8 iphase_not_imposed

        self.cP.iphase_critical_point
        Out: 4

        #if self.HEOS.phase() in [0,6]: # 0 is liquid phase
        #    self.HEOS.update(self.cP.QT_INPUTS, 1, TRef) # 1 = 100% vapor
        """
        self.HEOS = self.cP.AbstractState("HEOS", component)
        self.HEOS.specify_phase(self.cP.iphase_not_imposed)
        self._P = 101325
        self._T = 298.15
        self.cPFluids = dict([(cP.CoolProp.get_fluid_param_string(fluid,"CAS"),
                               fluid) for fluid in 
                              cP.CoolProp.get_global_param_string("fluids_list").split(',')])
        
    def vP(self,T):
        """
        vP(T)
        
        vapor pressure in Pa at T
        
        Parameters:
            T, temperature in K
        """
        self.HEOS.update(self.cP.QT_INPUTS, 1, T)
        return self.HEOS.p() 

    def vV(self,T):
        """
        vV(T)
        
        vapor viscosity in Pa*s at T
        
        Parameters:
            T, temperature in K
        """
        self.HEOS.update(self.cP.QT_INPUTS, 1, T)
        return self.HEOS.viscosity() 

    def lV(self,T):
        """
        lV(T)
        
        liquid viscosity in Pa*s at T
        
        Parameters:
            T, temperature in K
        """
        self.HEOS.update(self.cP.QT_INPUTS, 0, T)
        return self.HEOS.viscosity() 

    def lD(self,T):
        """
        lD(T)
        
        liquid mass density in kg/m^3 at T
        
        Parameters:
            T, temperature in K
        """
        self.HEOS.update(self.cP.QT_INPUTS, 0, T)
        return self.HEOS.rhomass() 

    def lTC(self,T):
        """
        lTC(T)
        
        liquid thermal conductivity in W/m/K at T
        
        Parameters:
            T, temperature in K
        """
        self.HEOS.update(self.cP.QT_INPUTS, 0, T)
        return self.HEOS.conductivity() 

    def vTC(self,T):
        """
        vTC(T)
        
        vapor thermal conductivity in W/m/K at T
        
        Parameters:
            T, temperature in K
        """
        self.HEOS.update(self.cP.QT_INPUTS, 1, T)
        return self.HEOS.conductivity() 

    def sT(self,T):
        """
        sT(T)
        
        liquid surface tension in N/m at T
        
        Parameters:
            T, temperature in K
        """
        self.HEOS.update(self.cP.QT_INPUTS, 0, T)
        return self.HEOS.surface_tension() 

    def eV(self,T):
        """
        eV(T)

        enthalpy of vaporization in kJ/mol at T
        
        Parameters:
            T, temperature in K

        """
        self.HEOS.update(self.cP.QT_INPUTS, 1, T)
        vH = self.HEOS.hmolar()
        self.HEOS.update(self.cP.QT_INPUTS, 0, T)
        lH = self.HEOS.hmolar()
        return (vH - lH)/1000 #kJ/mol
    
    def lHC(self,T):
        """
        lHC(T)

        liquid Heat Capacity in J/mol/K at T
        
        Parameters:
            T, temperature in K
            
        """        
        self.HEOS.update(self.cP.QT_INPUTS, 0, T)
        return self.HEOS.cpmolar()

    def lH(self,T,TRef=298.15):
        """
        lH(T,TRef=298.15)

        liquid enthalpy in J/mol at T
        
        Parameters:
            T, temperature in K
            TRef, temperature in K
            
        """        
        self.HEOS.update(self.cP.PT_INPUTS, P, T)
        #self.cP.CoolProp.get_phase_index('phase_liquid')
        if self.HEOS.phase() in [0,6]: # 0 is liquid, 5 is gas, 6 is 'twophase'
            self.HEOS.update(self.cP.QT_INPUTS, 0, T)
        lHT = self.HEOS.hmolar()
        self.HEOS.update(self.cP.QT_INPUTS, 0, TRef)
        lHTRef = self.HEOS.hmolar()
        return (lHT - lHTRef) #J/mol

    def vHC(self,T,P=101325):
        """
        vHC(T,P=101325)

        vapor Heat Capacity in J/mol/K at T 
        
        Parameters:
            T, temperature in K
            P, pressure in Pa
            
        """        
        self.HEOS.update(self.cP.PT_INPUTS, P, T)
        if self.HEOS.phase() in [0,6]: # 6 is 'twophase'
            self.HEOS.update(self.cP.QT_INPUTS, 1, T)
        return self.HEOS.cpmolar()

    def vH(self,T,P=1,TRef=298.15,PRef=1):
        """
        vH(T,P=1,TRef=298.15,PRef=1)

        vapor enthalpy in J/mol at T and P
        
        Parameters:
            T, temperature in K (T > freezing point at P = 1 Pa) 
            P, (partial) pressure of component in Pa
            TRef, temperature in K (T > freezing point at P = 1 Pa)
            PRef, (partial) pressure of component in Pa

        TODO:
            handle T< freezing point for low pressure for T only input.
            i.e. check/force T > db['water']['nFP']?
        """        
        self.HEOS.specify_phase(self.cP.iphase_gas)
        self.HEOS.update(self.cP.PT_INPUTS, P, T)
        H = self.HEOS.hmolar()
        self.HEOS.update(self.cP.PT_INPUTS, PRef, TRef)
        HRef = self.HEOS.hmolar()
        return (H - HRef) #J/mol

def enthalpy_f1(T, hCP, TRef=298.15):
    """
    enthalpy_f1(T, hCP, TRef=298.15)
      slop vapor & solid phases: T in K, enthalpy in J/mol;

    enthalpy change correlation
      enthalpy = A*(T-TRef) + 1/2*B*(T^2-TRef^2) - C*(1/T-1/TRef)
    
    Parameters
      T, temperature
      TRef, reference temperature
      hCP, A=hCP[0], B=hCP[1], C=hCP[2]
      A, B, and C are regression coefficients  

    Returns
      enthalpy change between T and TRef
    """
    return hCP[0]*(T-TRef) + 1.0/2.0*hCP[1]*(T**2-TRef**2) - hCP[2]*(1/T-1/TRef)

def enthalpy_Shomate(T, hCP):
    """
    enthalpy_Shomate(T, hCP)
      NIST vapor, liquid, and solid phases

    heat capacity correlation
      H - H_ref (kJ/mol) = A*t + 1/2*B*t^2 + 1/3*C*t^3 + 
                           1/4*D*t^4 - E*1/t + F - H
      t (K) = T/1000.0
      H_ref is enthalpy in kJ/mol at 298.15 K
    
    Parameters
      T, temperature in Kelvin
      hCP, A=hCP[0], B=hCP[1], C=hCP[2], D=hCP[3], E=hCP[4], F=hCP[5], H=hCP[7]
      A, B, C, D, E, F, and H are regression coefficients  

    Returns
      enthalpy in kJ/mol at T relative to enthalpy at 298.15 K
    """
    t = T/1000.0
    return hCP[0]*t + hCP[1]/2*t**2 + hCP[2]/3*t**3 + hCP[3]/4*t**4 - hCP[4]/t + hCP[5] - hCP[7]

def enthalpy_Shomate_bounded(T, TBounds, sPs, TRef=298.15):
    """
    enthalpy_Shomate_bounded(T, Tbounds, sPs, TRef=298.15)

    Shomate Equation for enthalpy using multiple temperature bounded 
    parameter sets.  The first parameter set is used for witch 
    Tmin <= T <= Tmax where Tmin and Tmax are the minimum and maximum 
    temperatures over which the parameter set is valid. 

    Parameters
      T, temperature in Kelvin
      TRef, reference temperature in Kelvin
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

   Returns
      enthalpy_Shomate(T, sP[index]) - HRef if T and TRef are within TBounds 
      (TRef = 298.15 K is always ok) otherwise None
    """
    HRef = 0
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

def heatCapacity_f1(T, hCP):
    """
    heatCapacity_form1(T, hCP)
      slop vapor & solid phases: T in K, heat capaity in J/mol/K;

    heat capacity correlation
      heat capacity = A + B*T + C/T^2
    
    Parameters
      T, temperature
      hCP, A=hCP[0], B=hCP[1], C=hCP[2]
      A, B, and C are regression coefficients  

    Returns
      heat capacity at T
    """
    return hCP[0] + hCP[1]*T + hCP[2]/T**2

def heatCapacity(T, hCP):
    """
    heatCapacity(T, hCP)

    heatCapacity = A + B*T + C*T^2
    
    Parameters
      T, temperature in Kelvin
      hCP, A=hCP[0], B=hCP[1], C=hCP[2]
      A, B, and C are regression coefficients  

    Returns
      heat capacity at T
    """
    return hCP[0] + hCP[1]*T + hCP[2]*T**2

def heatCapacity_Shomate(T, hCP):
    """
    heatCapacity_Shomate(T, hCP)
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
    t = T/1000.0
    return hCP[0] + hCP[1]*t + hCP[2]*t**2 + hCP[3]*t**3 + hCP[4]/t**2

def henry_rSSI(T, hbpSIP):
    """
    henery_rSSI(T, hbpSIP)

    Henry's law correlation in molality per pressure from R. Sander, SI units

    henery_rSSI (mol/(kg solvent)/Pa) = Ho*exp({dln(H)/d(1/T)}*(1/T-1/To))
      
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

def mCVTC(molFrac, vTC, vV, fW, eps=1.0):
    """
    mCVTC(molFrac, vTC, vV, fW, eps=1.0)
    
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
    tCm = 0
    Aij = []
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

def mCVV(molFrac, vV, fW):
    """
    mCVV(molFrac, vV, fW)
    
    multi-Component Vapor Viscosity in micropoise
    
    Parameters:
        molFrac, list of mol fractions
        vV, list of pure component vapor viscosities (components ordered as
            in fW mol)
        fW, list of formula weights for each component in molFrac
    
    Returns:
        multi-component vapor viscosity in micropoise
    """
    vm = 0
    Aij = []
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

def mFWQ(mFL, qL):
    """
    mFWQ(mFL, qL):
        
    (mole or mass) Fraction Weighted Quantity
    
    Parameters:
        mFL, list of mole or mass fractions, sum(mFL) = 1
        qL, list of quantities corresponding to items in mFL
    
    Returns:
        weighted averaged of items in qL
    """
    aveQ = 0
    for idx,mF in enumerate(mFL):
        aveQ += mF*qL[idx]
    
    return aveQ

def thermoFuncs(db={},useCoolProp=False):
    """
    add thermo functions from yaws_python.py, coolprop, etc to a db dictionary
    prefer CoolProp functions over yaws for all CoolProp fluids


    TODO: - CURRENTLY BROKEN while reorganizing "mpdb.eq"
          - fix mCM.py and other demos that depend on "mpdb.eq"
          - sensible strategy to use multiple EQs for same thermo functions
            i.e. CoolProp, pcmpEQ, mcmpEQ for vapor enthalpy...
          - coerce CoolProp into returning value for liquid or vapor states
            given only T like yaws eqns...
    NOTES:
          - python closures are late binding.  The value of variables used in 
            closures are looked up at the time an inner function is called (the 
            lambda's below)
          - python default function arguments are mutable and initialized once
            at the time the function is defined.  Mutalbe means the default 
            argument can be changed at each function call.

    """ 
    cPFluids = coolPropFluids()
    dbtf={}
    for k in db.keys():
        dbtf.update({k: {}})
        if useCoolProp:
            CAS = db[k]['CAS']
            if CAS in cPFluids.keys():
                cPK = cPFluids[CAS]
                dbtf[k].update({'cPKey': cPK})
                dbtf[k].update({'cPEQ': coolPropEqns(cPK)})
                
        if 'vPP' in db[k].keys():
            dbtf[k].update({'vP':lambda T, vPP=db[k]['vPP']:
                            yp.vaporPressure(T,vPP)})
        if 'eVP' in db[k].keys():
            dbtf[k].update({'eV':lambda T, eVP=db[k]['eVP']:
                            yp.enthalpyVaporization(T,eVP)})
        if 'lDP' in db[k].keys():
            dbtf[k].update({'lD':lambda T, lDP=db[k]['lDP']:
                            yp.liquidDensity(T,lDP)})
        if 'lHCP' in db[k].keys():
            dbtf[k].update({'lHC':lambda T, lHCP=db[k]['lHCP']:
                            yp.liquidHeatCapacity(T,lHCP)})
        if 'sHCP' in db[k].keys():
            dbtf[k].update({'sHC':lambda T, sHCP=db[k]['sHCP']:
                            yp.solidHeatCapacity(T,sHCP)})
            dbtf[k].update({'sH':lambda T, TRef, sHCP=db[k]['sHCP']:
                            yp.solidEnthalpy(T,sHCP,TRef=TRef)})
        if 'vHCP' in db[k].keys():
            dbtf[k].update({'vHC':lambda T, vHCP=db[k]['vHCP']:
                            yp.vaporHeatCapacity(T,vHCP)})
            dbtf[k].update({'vH':lambda T, TRef, vHCP=db[k]['vHCP']:
                            yp.vaporEnthalpy(T,vHCP,TRef=TRef)})
        """
        if 'hCP' in db[k].keys():
            dbtf[k].update({'sHC_f1':lambda T, sHCP=db[k]['sHCP_f1']: hC_f1(T,sHCP)})
            dbtf[k].update({'sH_f1':lambda T, TRef, sHCP=db[k]['sHCP_f1']: h_f1(T,sHCP,TRef=TRef)})
        if 'sHCP_f1p1' in db[k].keys():
            dbtf[k].update({'sHC_f1p1':lambda T, sHCP=db[k]['sHCP_f1p1']: hC_f1(T,sHCP)})
            dbtf[k].update({'sH_f1p1':lambda T, TRef, sHCP=db[k]['sHCP_f1p1']: h_f1(T,sHCP,TRef=TRef)})
        if 'sHCP_f1p2' in db[k].keys():
            dbtf[k].update({'sHC_f1p2':lambda T, sHCP=db[k]['sHCP_f1p2']: hC_f1(T,sHCP)})
            dbtf[k].update({'sH_f1p2':lambda T, TRef, sHCP=db[k]['sHCP_f1p2']: h_f1(T,sHCP,TRef=TRef)})
        if 'vHCP_f1' in db[k].keys():
            dbtf[k].update({'vHC_f1':lambda T, vHCP=db[k]['vHCP_f1']: hC_f1(T,vHCP)})
            dbtf[k].update({'vH_f1':lambda T, TRef, vHCP=db[k]['vHCP_f1']: h_f1(T,vHCP,TRef=TRef)})
        """
        if 'vVP' in db[k].keys():
            dbtf[k].update({'vV':lambda T, vVP=db[k]['vVP']:
                            yp.vaporViscosity(T,vVP)})
        if 'lVP' in db[k].keys():
            dbtf[k].update({'lV':lambda T, lVP=db[k]['lVP']:
                            yp.liquidViscosity(T,lVP)})
        if 'vTCP' in db[k].keys():
            dbtf[k].update({'vTC':lambda T, vTCP=db[k]['vTCP']:
                            yp.vaporThermalConductivity(T,vTCP)})
        if 'tCP' in db[k].keys():
            dbtf[k].update({'tC':lambda T, tCP=db[k]['tCP']:
                            yp.thermalConductivity(T,tCP)})
        if 'lTCP' in db[k].keys():
            dbtf[k].update({'lTC':lambda T, lTCP=db[k]['lTCP']:
                            yp.liquidThermalConductivity(T,lTCP)})
        if 'sTP' in db[k].keys():
            dbtf[k].update({'sT':lambda T, sTP=db[k]['sTP']:
                            yp.surfaceTension(T,sTP)})
        if 'vGEFP' in db[k].keys():
            dbtf[k].update({'vGEF':lambda T,vGEFP=db[k]['vGEFP']:
                            yp.vaporGibbsEnergyFormation(T,vGEFP)})
        if 'vEFP' in db[k].keys() and \
           'cPEQ' not in dbtf[k].keys():
            dbtf[k].update({'vEF':lambda T, vEFP=db[k]['vEFP']:
                            yp.vaporEnthalpyFormation(T,vEFP)})
        elif 'eFHUSState' in db[k].keys() and \
             'gas' in db[k]['eFHUSState'] and \
             'cPEQ' in dbtf[k].keys(): \
             dbtf[k].update({'vEF':lambda T, ck=k: db[ck]['eFH298K'] +\
                             dbtf[ck]['cPEQ'].vH(T)/1000})
        elif 'eFHUSState' in db[k].keys() and \
             'gas' in db[k]['eFHUSState'] and \
             'vHCP' in db[k].keys(): \
             dbtf[k].update({'vEF':lambda T, ck=k: db[ck]['eFH298K'] +\
                             dbtf[ck]['vH'](T,298.15)/1000})

    return dbtf
