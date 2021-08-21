"""
mpdb: Material Property Data Base as python module
misc_python.py: utility methods for thermophysical property estimation

Revision Date: 2021.08.17

SPDX-License-Identifier: BSD-2-Clause
Copyright (c) 2021 Stuart Nolan. All rights reserved.
"""
import CoolProp as cP

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
        Thermophysical equations similar to yaws equations in mpdb.pcmpEQ as
        calculated by CoolProp
        
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
        
        vapor pressure at T in Pa
        
        Parameters:
            T, temperature in K
        """
        self.HEOS.update(self.cP.QT_INPUTS, 1, T)
        return self.HEOS.p() 

    def vV(self,T):
        """
        vV(T)
        
        vapor viscosity at T in Pa*s
        
        Parameters:
            T, temperature in K
        """
        self.HEOS.update(self.cP.QT_INPUTS, 1, T)
        return self.HEOS.viscosity() 

    def lV(self,T):
        """
        lV(T)
        
        liquid viscosity at T in Pa*s
        
        Parameters:
            T, temperature in K
        """
        self.HEOS.update(self.cP.QT_INPUTS, 0, T)
        return self.HEOS.viscosity() 

    def lD(self,T):
        """
        lD(T)
        
        liquid mass density at T in kg/m^3
        
        Parameters:
            T, temperature in K
        """
        self.HEOS.update(self.cP.QT_INPUTS, 0, T)
        return self.HEOS.rhomass() 

    def lTC(self,T):
        """
        lTC(T)
        
        liquid thermal conductivity at T in W/m/K
        
        Parameters:
            T, temperature in K
        """
        self.HEOS.update(self.cP.QT_INPUTS, 0, T)
        return self.HEOS.conductivity() 

    def vTC(self,T):
        """
        vTC(T)
        
        vapor thermal conductivity at T in W/m/K
        
        Parameters:
            T, temperature in K
        """
        self.HEOS.update(self.cP.QT_INPUTS, 1, T)
        return self.HEOS.conductivity() 

    def sT(self,T):
        """
        sT(T)
        
        liquid surface tension at T in N/m
        
        Parameters:
            T, temperature in K
        """
        self.HEOS.update(self.cP.QT_INPUTS, 0, T)
        return self.HEOS.surface_tension() 

    def eV(self,T):
        """
        eV(T)

        enthalpy of vaporization at T in kJ/mol
        
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

        liquid Heat Capacity at T in J/mol/K
        
        Parameters:
            T, temperature in K
            
        """        
        self.HEOS.update(self.cP.QT_INPUTS, 0, T)
        return self.HEOS.cpmolar()

    def lH(self,T,TRef=298.15):
        """
        lH(T,TRef=298.15)

        liquid enthalpy at T in J/mol
        
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

        vapor Heat Capacity at T in J/mol/K
        
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

        vapor enthalpy at T and P in J/mol
        
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

def enthalpy_Shomate(T, hCP, TRef=298.15):
    """
    enthalpy_Shomate(T, hCP, TRef=298.15)
      NIST vapor, liquid, and solid phases

    heat capacity correlation
      hC (J/mol/K) = A*(t-tr) + 1/2*B*(t^2-tr^2) + 1/3*C*(t^3-tr^3) + 
                     1/4*D*(t^4-tr^4) - E*(1/t-1/tr)
      t (K) = T/1000.0
      tr (K) = TRef/1000.0
    
    Parameters
      T, temperature in Kelvin
      hCP, A=hCP[0], B=hCP[1], C=hCP[2], D=hCP[3], E=hCP[4]
      A, B, C, D, and E are regression coefficients  

    Returns
      heat capacity in J/mol/K at T
    """
    t = T/1000.0
    tr = TRef/1000.0
    return hCP[0]*(t-tr) + hCP[1]*(t**2-tr**2)/2 + hCP[2]*(t**3-tr**3)/3 +
           hCP[3]*(t**4-t**4)/4 - hCP[4]*(1/t-1/tr)

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


