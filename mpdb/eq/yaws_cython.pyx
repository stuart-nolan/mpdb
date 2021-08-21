"""
mpdb: Material Property Data Base as python module
yaws_cython.pyx: equations for "yaws" thermophysical proptery estimation

Revision Date: 2021.08.16

SPDX-License-Identifier: BSD-2-Clause
Copyright (c) 2021 Stuart Nolan. All rights reserved.

To build:
python yaws_cython-setup.py build_ext --inplace
"""

from libc.math cimport log10, exp

cimport cython
@cython.boundscheck(False) # turn off bounds-checking
@cython.wraparound(False)  # turn off negative index wrapping

cpdef float activatedCarbonAdsorption(float y, list aCAP):
    """
    activatedCarbonAdsorption(float y, list aCAP)

    activatedCarbonAdsorption (g/100 g of carbon) = 10**(A + B*log10(y) + 
                                                    C*log10(y)**2)
    
    Parameters
      y, concentration of gas at 25 Deg. C at 1 atm in ppm by volume
      aCAP, A=aCAP[0], B=aCAP[1], C=aCAP[2]
      A, B, and C are regression coefficients  

    Returns
      adsorption capactiy at equilibirum in g/100 g carbon
    """
    cdef float log10y
    
    log10y = log10(y)
    return 10**(aCAP[0] + aCAP[1]*log10y + aCAP[2]*log10y**2)

cpdef float enthalpyVaporization(float T, list eVP):
    """
    enthalpyVaporization(float T, list eVP)
    
    enthalpyVaporization (kJ/mol) = A*(1-T/critT)^n
    
    Parameters
      T, temperature in Kelvin
      eVP, A=eVP[0], critT=eVP[1], n=eVP[2]
      A and n are regression coefficients, critT: critical temperature

    Returns
      enthalpy of vaporization in kJ/mol at T
    """
    return eVP[0]*(1-T/eVP[1])**eVP[2]
        
cpdef float liquidDensity(float T, list lDP):
    """
    liquidDensity(float T, list lDP)
    
    liquidDensity (kg/m^3) = A*B^-(1-T/critT)^n*1000.0
    
    Parameters
      T, temperature in Kelvin
      lDP, A=lDP[0], B=lDP[1], n=lDP[2], critT=lDP[3]
      A, B, and n are regression coefficients, critT: critical temperature

    Returns
      saturated liquid density in kg/m^3 at T
    """
    return 1000.0*lDP[0]*lDP[1]**-((1-T/lDP[3])**lDP[2])

cpdef float liquidHeatCapacity(float T, list lHCP):
    """
    liquidHeatCapacity(float T, list lHCP)
    
    liquidHeatCapacity (J/mol/K) = A + B*T + C*T^2 + D*T^3
    
    Parameters
      T, temperature in Kelvin
      lHCP, A=lHCP[0], B=lHCP[1], C=lHCP[2], D=lHCP[3]
      A, B, C, and D are regression coefficients  

    Returns
      liquid heat capacity in J/mol/K at T
    """
    return lHCP[0] + lHCP[1]*T + lHCP[2]*T**2 + lHCP[3]*T**3

cpdef float liquidEnthalpy(float T, list lHCP, float TRef=298.15):
    """
    liquidEnthalpy(float T, list lHCP, float TRef=298.15)
    
    liquidEnthalpy (J/mol) = A*(T-TRef) + 1/2*B*(T^2-TRef^2) + 
                             1/3*C*(T^3-TRef^3) + 1/4*D*(T^4-TRef^4)
    
    Parameters
      T, temperature in Kelvin
      TRef, reference temperature in Kelvin
      lHCP, A=lHCP[0], B=lHCP[1], C=lHCP[2], D=lHCP[3]
      A, B, C, and D are regression coefficients  

    Returns
      liquid enthalpy change in J/mol between T and TRef
    """
    return lHCP[0]*(T-TRef) + 1.0/2.0*lHCP[1]*(T**2-TRef**2) + 1.0/3.0*lHCP[2]*(T**3-TRef**3) + 1.0/4.0*lHCP[3]*(T**4-TRef**4)

cpdef float liquidThermalConductivity(float T, list lTCP):
    """
    liquidThermalConductivity(float T, list lTCP)

    liquidThermalConductivity (W/m/K) = 10**(A + B*(1-T/C)^(2/7))
    
    Parameters
      T, temperature in Kelvin
      lTCP, A=lTCP[0], B=lTCP[1], C=lTCP[2]
      A, B, and C are regression coefficients  

    Returns
      (organic) liquid thermal conductivity in W/m/K at T
    """
    return 10**(lTCP[0] + lTCP[1]*(1-T/lTCP[2])**(2.0/7.0))

cpdef float liquidViscosity(float T, list lVP):
    """
    liquidViscosity(float T, list lVP)
    
    liquidViscosity (centipoise) = 10^(A + B/T + C*T + D*T^2)
    
    Parameters
      T, temperature in K
      vPP, A=lVP[0], B=lVP[1], C=lVP[2], D=lVP[3]
      A, B, C, D and E are regression coefficients
      
    Returns
      liquid viscosity in centipoise at T
    """  
    return 10**(lVP[0] + lVP[1]/T + lVP[2]*T + lVP[3]*T**2)

cpdef float saltWaterSolubility(float X, list sWSP):
    """
    saltWaterSolubility(float X, list sWSP)

    saltWaterSolubility (ppm by wt) = 10**(A + B*X + C*X**2)
    
    Parameters
      X, concentration of NaCl in water at 25 Deg. C in ppm by wt
      sWSP, A=sWSP[0], B=sWSP[1], C=sWSP[2]
      A, B, and C are regression coefficients  

    Returns
      salt water solubility in ppm by weight at X
    """
    return 10**(sWSP[0] + sWSP[1]*X + sWSP[2]*X**2)

cpdef float solidHeatCapacity(float T, list sHCP):
    """
    solidHeatCapacity(float T, list sHCP)

    solidHeatCapacity (J/mol/K) = A + B*T + C*T^2
    
    Parameters
      T, temperature in Kelvin
      sHCP, A=sHCP[0], B=sHCP[1], C=sHCP[2]
      A, B, and C are regression coefficients  

    Returns
      solid heat capacity in J/mol/K at T
    """
    return sHCP[0] + sHCP[1]*T + sHCP[2]*T**2

cpdef float solidEnthalpy(float T, list sHCP, float TRef=298.15):
    """
    solidEnthalpy(float T, list sHCP, float TRef=298.15)

    solidEnthalpy (J/mol) = A*(T-TRef) + 1/2*B*(T^2-TRef^2) + 
                            1/3*C*(T^3-TRef^3)
    
    Parameters
      T, temperature in Kelvin
      TRef, reference temperature in Kelvin
      sHCP, See solidHeatCapacity: A=sHCP[0], B=sHCP[1], C=sHCP[2]
      A, B, and C are regression coefficients  

    Returns
      solid enthalpy change in J/mol between T and TRef
    """
    return sHCP[0]*(T-TRef) + 1.0/2.0*sHCP[1]*(T**2-TRef**2) + 1.0/3.0*sHCP[2]*(T**3-TRef**3)

cpdef float surfaceTension(float T, list sTP):
    """
    surfaceTension(float T, list sTP)
    
    surfaceTension (dyne/cm) = A*(1-T/critT)^n
    
    Parameters
      T, temperature in Kelvin
      sTP, A=sTP[0], critT=sTP[1], n=sTP[2]
      A and n are regression coefficients, critT: critical temperature

    Returns
      surface tension in dyne/cm at T
    """
    return sTP[0]*(1-T/sTP[1])**sTP[2]

cpdef float thermalConductivity(float T, list tCP):
    """
    thermalConductivity(float T, list tCP)

    thermalConductivity (W/m/K) = A + B*T + C*T^2
    
    Parameters
      T, temperature in Kelvin
      tCP, A=tCP[0], B=tCP[1], C=tCP[2]
      A, B, and C are regression coefficients  

    Returns
      inorganic liquid or solid thermal conductivity in W/m/K at T
    """
    return tCP[0] + tCP[1]*T + tCP[2]*T**2

cpdef float vaporEnthalpy(float T, list vHCP, float TRef=298.15):
    """
    vaporEnthalpy(float T, list vHCP, float TRef=298.15)

    vaporEnthalpy (J/mol) = A*(T-TRef) + 1/2*B*(T^2-TRef^2) + 
                            1/3*C*(T^3-TRef^3) + 1/4*D*(T^4-TRef^4) + 
                            1/5*E*(T^5-TRef^5)

    Parameters
      T, temperature in Kelvin
      TRef, reference temperature in Kelvin
      vHCP, A=vHCP[0], B=vHCP[1], C=vHCP[2], D=vHCP[3], E=vHCP[4]
      A, B, C, D and E are regression coefficients  

    Returns
      vapor enthalpy change in J/mol between T and TRef
    """
    return vHCP[0]*(T-TRef) + 1.0/2.0*vHCP[1]*(T**2-TRef**2) + 1.0/3.0*vHCP[2]*(T**3-TRef**3) + 1.0/4.0*vHCP[3]*(T**4-TRef**4) + 1.0/5.0*vHCP[4]*(T**5-TRef**5)

cpdef float vaporEnthalpyFormation(float T, list vEFP):
    """
    vaporEnthalpyFormation(float T, list vEFP)

    vaporEnthalpyFormation (kJ/mol) = A + B*T + C*T^2
    
    Parameters
      T, temperature in Kelvin
      vEFP, A=vEFP[0], B=vEFP[1], C=vEFP[2]
      A, B, and C are regression coefficients  

    Returns
      vapor Enthalpy of Formation in kJ/mol at T
    """
    return vEFP[0] + vEFP[1]*T + vEFP[2]*T**2

cpdef float vaporGibbsEnergyFormation(float T, list vGEFP):
    """
    vaporGibbsEnergyFormation(float T, list vGEFP)

    vaporGibbsEnergyFormation (kJ/mol) = A + B*T + C*T^2
    
    Parameters
      T, temperature in Kelvin
      vGEFP, A=vGEFP[0], B=vGEFP[1], C=vGEFP[2]
      A, B, and C are regression coefficients  

    Returns
      vapor Gibbs Energy of Formation in kJ/mol at T
    """
    return vGEFP[0] + vGEFP[1]*T + vGEFP[2]*T**2

cpdef float vaporHeatCapacity(float T, list vHCP):
    """
    vaporHeatCapacity(float T, list vHCP)

    vaporHeatCapacity (J/mol/K) = A + B*T + C*T^2 + D*T^3 + E*T^4

    Parameters
      T, temperature in Kelvin
      vHCP, A=vHCP[0], B=vHCP[1], C=vHCP[2], D=vHCP[3], E=vHCP[4]
      A, B, C, D and E are regression coefficients  

    Returns
      vapor heat capacity in J/mol/K at T
    """
    return vHCP[0] + vHCP[1]*T + vHCP[2]*T**2 + vHCP[3]*T**3 + vHCP[4]*T**4

cpdef float vaporPressure(float T, list vPP):
    """
    vaporPressure(float T, list vPP)
    
    vaporPressure (mmHg) = 10^(A + B/T + C*log10(T) + D*T + E*T^2)
    
    Parameters
      T, temperature in K
      vPP, A=vPP[0], B=vPP[1], C=vPP[2], D=vPP[3], E=vPP[4]
      A, B, C, D and E are regression coefficients  
      
    Returns
      vapor pressure in mmHg at T
    """  
    return 10**(vPP[0] + vPP[1]/T + vPP[2]*log10(T) + vPP[3]*T + vPP[4]*T**2)
            
cpdef float vaporThermalConductivity(float T, list vTCP):
    """
    vaporThermalConductivity(float T, list vTCP)

    vaporThermalConductivity (W/m/K) = A + B*T + C*T^2
    
    Parameters
      T, temperature in Kelvin
      vTCP, A=vTCP[0], B=vTCP[1], C=vTCP[2]
      A, B, and C are regression coefficients  

    Returns
      vapor thermal conductivity in W/m/K at T
    """
    return vTCP[0] + vTCP[1]*T + vTCP[2]*T**2

cpdef float vaporViscosity(float T, list vVP):
    """
    vaporViscosity(float T, list vVP)

    vaporViscosity (micropoise) = A + B*T + C*T^2
    
    Parameters
      T, temperature in Kelvin
      vVP, A=vVP[0], B=vVP[1], C=vVP[2]
      A, B, and C are regression coefficients  

    Returns
      vapor viscosity in micropoise at T
    """
    return vVP[0] + vVP[1]*T + vVP[2]*T**2

cpdef float waterSolubility(float T, list wSP):
    """
    waterSolubility(float T, list wSP)

    waterSolubility (ppm by wt) = 10**(A + B/T + C/T^2)
    
    Parameters
      T, temperature in Kelvin
      wSP, A=wSP[0], B=wSP[1], C=wSP[2]
      A, B, and C are regression coefficients  

    Returns
      water solubility in ppm by weight at T
    """
    return 10**(wSP[0] + wSP[1]/T + wSP[2]/T**2)

def yawsEqn(obj,method,*args):
    """
    yaws equation cython function dispatcher

    usage:

    from mpdb.eq import yawsEqns as yeqs
    from mpdb.utils import dbLoad
    (db, dbmd) = dbLoad(interactive=True)
    vP = yeqs.yawsEqn(yeqs,"vaporPressure")
    #vP = yeqs.vaporPressure # alternatively
    vP(298.15,db['water']['vPP'])
    23.782012939453125
    yeqs.vaporPressure(298.15,db['water']['vPP'])
    23.782012939453125
    %timeit yeqs.vaporPressure(298.15,db['water']['vPP'])
    523 ns ± 1.57 ns per loop (mean ± std. dev. of 7 runs, 1000000 loops each)
    %timeit vP(298.15,db['water']['vPP'])
    526 ns ± 7.11 ns per loop (mean ± std. dev. of 7 runs, 1000000 loops each)
    """
    return getattr(obj,'%s' % method)
