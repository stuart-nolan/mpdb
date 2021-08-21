"""
mpdb: Material Property Data Base as python package
units.py: wrapper functions to nix/gnuwin32 "units" CLI program

Revision Date: 2021.08.15

SPDX-License-Identifier: BSD-2-Clause
Copyright (c) 2021 Stuart Nolan. All rights reserved.

ALTERNATIVES:
    https://github.com/TheGrum/python3-physics
    https://bitbucket.org/birkenfeld/ipython-physics
    https://pint.readthedocs.io/en/latest/
    https://pythonhosted.org/uncertainties/
    https://pypi.org/project/quantities/
    https://pypi.org/project/units/
    http://docs.enthought.com/scimath/units/unit_numpy.html
"""
import os
    
def uc(value, inUnit, outUnit):
    """
    uc(value, inUnit, outUnit)
    
    Wrapper for the units(.exe) command.  Converts value in inUnits
    to outUnits.

    Parameters
    - value, numerical value in inUnits
    - inUnit, string representation of units to convert from
    - outUnit, string representation of units to convert to

    Returns
    - conversion from units command

    Examples
    - print("55 mph is %s km/hour" % uc(55, 'mph', 'km/hour'))
    - print("gas constant: %s J/mol/K" % uc(1,"R","J/mol/K"))
    - print("Boltzman's constant: %s J/K" % uc(1,"k","J/K"))
    - print("Avogadro's number: %s molecules/mol" % uc(1,"N_A","1/mol"))
    - print("pi: %s " % uc(1,"pi",""))
    - print("Faraday: %s C(/mol)" % uc(1,"faraday","C"))
    - print("electron charge: %s C" % uc(1,"e","C"))
    - print("A 180 deg. F temp. difference is a %s deg. C temp. difference" 
            % uc(180,'degF','degC'))

    """
    cmd = "/usr/bin/units -t"

    if not isinstance(value,(int,float)):
        return "uc: %s is an invalid value" % value
    if value < 0:
        value = "%s*-1" % (-value) #units command interprets "-" as an option

    try:
        res = os.popen("%s \"%s %s\" \"%s\"" % (cmd, value, inUnit, outUnit)).read()
        return float(res)
    except:
        return "uc: failed with: %s" % res
        
def ucTemp(value, inUnit, outUnit):
    """
    ucTemp(value, inUnit, outUnit)
    
    Wrapper for the units(.exe) command specific to absolute temperature 
    conversions.  

    Parameters
    - value, numerical value in inUnits
    - inUnit, string representation of temperature units to convert from
    - outUnit, string representation of temperature units units to convert to

    Valid entries for inUnit and outUnit are C, F, K, R (case sensitive)

    Returns
    - absolute temperature conversion from units command

    Examples
    - print("140 deg. F is %s deg. C" % ucTemp(140,'F','C'))

    """
    cmd = "/usr/bin/units -t"
    validT = ('C','K','F','R')
    if inUnit not in validT:
        return "ucTemp: %s is an invalid inUnit" % inUnit
    if outUnit not in validT:
        return "ucTemp: %s is an invalid outUnit" % outUnit
    if not isinstance(value,(int,float)):
        return "ucTemp: %s is an invalid value" % value
    inUnit = "temp%s(%s)" % (inUnit,value)
    outUnit = "temp" + outUnit
    try:
        res = os.popen("%s \"%s\" \"%s\"" % (cmd, inUnit, outUnit)).read()
        return float(res)
    except:
        return "ucTemp failed with: %s" % res

class gasConstant:
    def __init__(self,unit="J/mol/K"):
        """
        gas constant
        
        Parameters:
            unit, string like "cal/mol/degC"

        Usage:
            gC = gasConstant("psi*ft^3/mol/degF")
            print("gas constant in %s: %s" % (gC.unit,gC.R))
            gC.unit = "Pa*m^3/mol/K"
            print("gas constant in %s: %s" % (gC.unit,gC.R))
        """
        self._unit = "J/mol/K"
        self.R = uc(1,"R",self._unit)
        self.unit = unit

    @property
    def unit(self):
        return self._unit

    @unit.setter
    def unit(self,unit):
        R = uc(1,"R",unit)
        if type(R) is float:
            self._unit = unit
            self.R = R
        else:
            print("gasConstant: unit (%s) ERROR: %s; aborted" % (unit,R))
        
if __name__ == "__main__":
    print("55 mph is %s km/hour" % uc(55, 'mph', 'km/hour'))
    print("gas constant: %s J/mol/K" % uc(1,"R","J/mol/K"))
    print("Boltzman's constant: %s J/K" % uc(1,"k","J/K"))
    print("Avogadro's number: %s molecules/mol" % uc(1,"N_A","1/mol"))
    print("pi: %s " % uc(1,"pi",""))
    print("Faraday: %s C(/mol)" % uc(1,"faraday","C"))
    print("electron charge: %s C" % uc(1,"e","C"))
    print("A 180 deg. F temp. difference is a %s deg. C temp. difference" % uc(180,'degF','degC'))
    print("140 deg. F is %s deg. C" % ucTemp(140,'F','C'))

    gC = gasConstant("psi*ft^3/mol/degF")
    print("gas constant in %s: %s" % (gC.unit,gC.R))
    gC.unit = "Pa*m^3/mol/K"
    print("gas constant in %s: %s" % (gC.unit,gC.R))
