"""
mpdb: Material Property Data Base as python module
yaws_python.py: yaws equations for thermophysical property estimation

Revision Date: 2021.08.17

SPDX-License-Identifier: BSD-2-Clause
Copyright (c) 2021 Stuart Nolan. All rights reserved.
"""
from math import log10, exp

def activatedCarbonAdsorption(y, aCAP):
    """
    activatedCarbonAdsorption(y, aCAP)

    activatedCarbonAdsorption (g/100 g of carbon) = 10**(A + B*log10(y) + 
                                                    C*log10(y)**2)
    
    Parameters
      y, concentration of gas at 25 Deg. C at 1 atm in ppm by volume
      aCAP, A=aCAP[0], B=aCAP[1], C=aCAP[2]
      A, B, and C are regression coefficients  

    Returns
      adsorption capactiy at equilibirum in g/100 g carbon
    """
    log10y = log10(y)
    return 10**(aCAP[0] + aCAP[1]*log10y + aCAP[2]*log10y**2)

def enthalpyVaporization(T, eVP):
    """
    enthalpyVaporization(T, eVP)
    
    enthalpyVaporization (kJ/mol) = A*(1-T/critT)^n
    
    Parameters
      T, temperature in Kelvin
      eVP, A=eVP[0], critT=eVP[1], n=eVP[2]
      A and n are regression coefficients, critT: critical temperature

    Returns
      enthalpy of vaporization in kJ/mol at T
    """
    return eVP[0]*(1-T/eVP[1])**eVP[2]
        
def liquidDensity(T, lDP):
    """
    liquidDensity(T, lDP)
    
    liquidDensity (kg/m^3) = A*B^-(1-T/critT)^n*1000.0
    
    Parameters
      T, temperature in Kelvin
      lDP, A=lDP[0], B=lDP[1], n=lDP[2], critT=lDP[3]
      A, B, and n are regression coefficients, critT: critical temperature

    Returns
      liquid density in kg/m^3 at T
    """
    return 1000.0*lDP[0]*lDP[1]**-((1-T/lDP[3])**lDP[2])

def liquidHeatCapacity(T, hCP):
    """
    liquidHeatCapacity(T, hCP)

    liquidHeatCapacity = A + B*T + C*T^2 + D*T^3
    
    Parameters
      T, temperature
      hCP, A=hCP[0], B=hCP[1], C=hCP[2], D=hCP[3]
      A, B, C, and D are regression coefficients  

    Returns
      heat capacity at T
    """
    return lHCP[0] + lHCP[1]*T + lHCP[2]*T**2 + lHCP[3]*T**3

def liquidEnthalpy(T, hCP, TRef=298.15):
    """
    liquidEnthalpy(T, hCP, TRef=298.15)

    liquidEnthalpy (J/mol) = A*(T-TRef) + 1/2*B*(T^2-TRef^2) + 
                             1/3*C*(T^3-TRef^3) + 1/4*D*(T^4-TRef^4)
    
    Parameters
      T, temperature
      TRef, reference temperature
      hCP, A=hCP[0], B=hCP[1], C=hCP[2], D=hCP[3]
      A, B, C, and D are regression coefficients  

    Returns
      liquid enthalpy change in J/mol between T and TRef
    """
    return hCP[0]*(T-TRef) + 1.0/2.0*hCP[1]*(T**2-TRef**2) + 1.0/3.0*hCP[2]*(T**3-TRef**3) + 1.0/4.0*hCP[3]*(T**4-TRef**4)

def liquidThermalConductivity(T, lTCP):
    """
    liquidThermalConductivity(float T, list lTCP)

    liquidThermalConductivity (W/m/K) = 10**(A + B*(1-T/C)^(2/7))
    
    Parameters
      T, temperature in Kelvin
      lTCP, A=lTCP[0], B=lTCP[1], C=lTCP[2]
      A, B, and C are regression coefficients  

    Returns
      liquid thermal conductivity in W/m/K at T
    """
    return 10**(oLTCP[0] + oLTCP[1]*(1-T/oLTCP[2])**(2.0/7.0))

def liquidViscosity(T, lVP):
    """
    liquidViscosity(T, lVP)
    
    liquidViscosity (centipoise) = 10^(A + B/T + C*T + D*T^2)
    
    Parameters
      T, temperature in K
      vPP, A=lVP[0], B=lVP[1], C=lVP[2], D=lVP[3]
      A, B, C, D and E are regression coefficients
      
    Returns
      liquid viscosity in centipoise at T
    """  
    return 10**(lVP[0] + lVP[1]/T + lVP[2]*T + lVP[3]*T**2)

def saltWaterSolubility(X, sWSP):
    """
    saltWaterSolubility(X, sWSP)

    saltWaterSolubility (ppm by wt) = 10**(A + B*X + C*X**2)
    
    Parameters
      X, concentration of NaCl in water at 25 Deg. C in ppm by wt
      sWSP, A=sWSP[0], B=sWSP[1], C=sWSP[2]
      A, B, and C are regression coefficients  

    Returns
      salt water solubility in ppm by weight at X
    """
    return 10**(sWSP[0] + sWSP[1]*X + sWSP[2]*X**2)

def solidHeatCapacity(T, sHCP):
    """
    solidHeatCapacity(T, sHCP)

    solidHeatCapacity (J/mol/K) = A + B*T + C*T^2
    
    Parameters
      T, temperature in Kelvin
      sHCP, A=sHCP[0], B=sHCP[1], C=sHCP[2]
      A, B, and C are regression coefficients  

    Returns
      solid heat capacity in J/mol/K at T
    """
    return sHCP[0] + sHCP[1]*T + sHCP[2]*T**2

def solidEnthalpy(T, hCP, TRef=298.15):
    """
    solidEnthalpy(T, hCP, TRef=298.15)

    solidEnthalpy (J/mol) = A*(T-TRef) + 1/2*B*(T^2-TRef^2) + 
                            1/3*C*(T^3-TRef^3)
    
    Parameters
      T, temperature in Kelvin
      TRef, reference temperature
      hCP, A=hCP[0], B=hCP[1], C=hCP[2]
      A, B, and C are regression coefficients  

    Returns
      enthalpy change in J/mol between T and TRef
    """
    return hCP[0]*(T-TRef) + 1.0/2.0*hCP[1]*(T**2-TRef**2) + 1.0/3.0*hCP[2]*(T**3-TRef**3)

def surfaceTension(T, sTP):
    """
    surfaceTension(T, sTP)
    
    surfaceTension (dyne/cm) = A*(1-T/critT)^n
    
    Parameters
      T, temperature in Kelvin
      sTP, A=sTP[0], critT=sTP[1], n=sTP[2]
      A and n are regression coefficients, critT: critical temperature

    Returns
      surface tension in dyne/cm at T
    """
    return sTP[0]*(1-T/sTP[1])**sTP[2]

def thermalConductivity(T, tCP):
    """
    thermalConductivity(T, tCP)

    thermalConductivity (W/m/K) = A + B*T + C*T^2
    
    Parameters
      T, temperature in Kelvin
      tCP, A=tCP[0], B=tCP[1], C=tCP[2]
      A, B, and C are regression coefficients  

    Returns
      thermal conductivity in W/m/K at T
    """
    return iOTCP[0] + iOTCP[1]*T + iOTCP[2]*T**2

def vaporEnthalpy(T, hCP, TRef=298.15):
    """
    vaporEnthalpy(T, hCP, TRef=298.15)

    vaporEnthalpy (J/mol) = A*(T-TRef) + 1/2*B*(T^2-TRef^2) + 
                            1/3*C*(T^3-TRef^3) + 1/4*D*(T^4-TRef^4) + 
                            1/5*E*(T^5-TRef^5)

    Parameters
      T, temperature in K
      TRef, reference temperature
      hCP, A=hCP[0], B=hCP[1], C=hCP[2], D=hCP[3], E=hCP[4]
      A, B, C, D and E are regression coefficients  

    Returns
      enthalpy change in J/mol between T and TRef
    """
    return hCP[0]*(T-TRef) + 1.0/2.0*hCP[1]*(T**2-TRef**2) + 1.0/3.0*hCP[2]*(T**3-TRef**3) + 1.0/4.0*hCP[3]*(T**4-TRef**4) + 1.0/5.0*hCP[4]*(T**5-TRef**5)

def vaporEnthalpyFormation(T, vEFP):
    """
    vaporEnthalpyFormation(T, vEFP)

    vaporEnthalpyFormation (kJ/mol) = A + B*T + C*T^2
    
    Parameters
      T, temperature in Kelvin
      vEFP, A=vEFP[0], B=vEFP[1], C=vEFP[2]
      A, B, and C are regression coefficients  

    Returns
      vapor Enthalpy of Formation in kJ/mol at T
    """
    return vEFP[0] + vEFP[1]*T + vEFP[2]*T**2

def vaporHeatCapacity(T, hCP):
    """
    vaporHeatCapacity(T, hCP)

    vaporHeatCapacity (J/mol) = A + B*T + C*T^2 + D*T^3 + E*T^4

    Parameters
      T, temperature in K
      hCP, A=hCP[0], B=hCP[1], C=hCP[2], D=hCP[3], E=hCP[4]
      A, B, C, D and E are regression coefficients

    Returns
      heat capacity in J/mol at T
    """
    return hCP[0] + hCP[1]*T + hCP[2]*T**2 + hCP[3]*T**3 + hCP[4]*T**4

def vaporGibbsEnergyFormation(T, vGEFP):
    """
    vaporGibbsEnergyFormation(T, vGEFP)

    vaporGibbsEnergyFormation (kJ/mol) = A + B*T + C*T^2
    
    Parameters
      T, temperature in Kelvin
      vGEFP, A=vGEFP[0], B=vGEFP[1], C=vGEFP[2]
      A, B, and C are regression coefficients  

    Returns
      vapor Gibbs Energy of Formation in kJ/mol at T
    """
    return vGEFP[0] + vGEFP[1]*T + vGEFP[2]*T**2

def vaporPressure(T, vPP):
    """
    vaporPressure(T, vPP)
    
    vaporPressure (mmHg) = 10^(A + B/T + C*log10(T) + D*T + E*T^2)
    
    Parameters
      T, temperature in K
      vPP, A=vPP[0], B=vPP[1], C=vPP[2], D=vPP[3], E=vPP[4]
      A, B, C, D and E are regression coefficients  
      
    Returns
      vapor pressure in mmHg at T
    """  
    return 10**(vPP[0] + vPP[1]/T + vPP[2]*log10(T) + vPP[3]*T + vPP[4]*T**2)
            
def vaporThermalConductivity(T, vTCP):
    """
    vaporThermalConductivity(T, vTCP)

    vaporThermalConductivity (W/m/K) = A + B*T + C*T^2
    
    Parameters
      T, temperature in Kelvin
      vTCP, A=vTCP[0], B=vTCP[1], C=vTCP[2]
      A, B, and C are regression coefficients  

    Returns
      vapor thermal conductivity in W/m/K at T
    """
    return vTCP[0] + vTCP[1]*T + vTCP[2]*T**2

def vaporViscosity(T, vP):
    """
    vaporViscosity(T, vP)

    vaporViscosity (micropoise) = A + B*T + C*T^2
    
    Parameters
      T, temperature in Kelvin
      vVP, A=vVP[0], B=vVP[1], C=vVP[2]
      A, B, and C are regression coefficients  

    Returns
      vapor viscosity in micropoise at T
    """
    return vVP[0] + vVP[1]*T + vVP[2]*T**2

def waterSolubility(T, wSP):
    """
    waterSolubility(T, wSP)

    waterSolubility (ppm by wt) = 10**(A + B/T + C/T^2)
    
    Parameters
      T, temperature in Kelvin
      wSP, A=wSP[0], B=wSP[1], C=wSP[2]
      A, B, and C are regression coefficients  

    Returns
      water solubility in ppm by weight at T
    """
    return 10**(wSP[0] + wSP[1]/T + wSP[2]/T**2)
