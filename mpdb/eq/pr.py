"""
mpdb: python Material Property Data Base
pr.py: Ping Robinson EOS.  Partial implmentation mostly abandoned in favor of
the EOS availabe in CoolProp.

Revision Date: 2021.08.15

SPDX-License-Identifier: BSD-2-Clause
Copyright (c) 2021 Stuart Nolan. All rights reserved.
"""
from math import log
from math import exp
from math import sqrt
from mpdb.Constants import *
from mpdb.Units import *
from scipy import *

class PCPR:
    def __init__(self, cmpd, method="Pressure"):
        cst = Const.Const()
        self.R = cst.R["SI"]
        self.PUnit = "Pa"
        self.TUnit = "K"
        self.vUnit = "m^3/g-mol"
        self.Pc = units.uc(cmpd["CritP"]["Value"],\
                           cmpd["CritP"]["VUnit"],\
                           self.PUnit)
        self.Tc = units.uc(cmpd["CritT"]["Value"],\
                           cmpd["CritT"]["VUnit"],\
                           self.TUnit)
        self.omega = cmpd["Omega"]
        self.b = 0.07780*self.R*self.Tc/self.Pc
        self.fomega = 0.37464 + 1.54226*self.omega - 0.26992*self.omega**2
        
        self.method = method
        
    def value(self, *args):
        """
        """
        PCPR_Method = getattr(self, "PCPR_%s" % self.method)

        return PCPR_Method(*args)

    def PCPR_Pressure(self, T, v):
        """Pressure from Peng-Robinson equation of state.

        Parameters
          - T, temperature in kelvin
          - v, molar volume in m^3/g-mol

        Returns
          - P = RT/(v-b) - a/(v^2 + 2bv - b^2)
        where
           b = 0.07780*R*Tc/Pc
           a = 0.45724*R^2*Tc^2*[1+fw*(1-Tr^0.5)]^2/Pc
           fw = 0.337464 + 1.54226*w - 0.26992*w^2

           """
        a = 0.45724*self.R**2*self.Tc**2* \
            (1+self.fomega*(1-(T/self.Tc)**0.5))**2/self.Pc
        return self.R*T/(v-self.b) - a/(v**2 + 2*self.b*v - self.b**2)

    def PCPR_FugacityCoef(self, T, P, Z):
        """Fugacity coefficient from Peng-Robinson equation of state.

        Parameters
          - T, temperature in kelvin
          - P, pressure in Pa
          - Z, compressibility

        Returns
          - fugacity in Pa
          """

        a = 0.45724*self.R**2*self.Tc**2* \
            (1+self.fomega*(1-(T/self.Tc)**0.5))**2/self.Pc
        PbRT = P*self.b/self.R/T
        sqrt2 = sqrt(2)
        return exp(Z - 1 - log(Z-PbRT) - a*log((Z+(1+sqrt2)*PbRT)/ \
                                                  (Z+(1-sqrt2)*PbRT))/ \
                   2/sqrt2/self.b/self.R/T)

    def PCPR_Z(self, T, P):
        """Returns the liquid and vapor phase compressibility's

        Parameters
          - P, pressure in Pa
          - T, temperature in kelvin

        Returns
          - list of gas and/or liquid compresibilities
            - for VLE: 
              [ZV, ZL]
            - for single phase liquid or gas
              [Z]
          """
        a = 0.45724*self.R**2*self.Tc**2* \
            (1+self.fomega*(1-(T/self.Tc)**0.5))**2/self.Pc
        A = a*P/self.R**2/T**2
        B = P*self.b/self.R/T
        a0 = B**3+B**2-A*B
        a1 = A-3*B**2-2*B
        a2 = B-1
        p = (3*a1-a2**2)/3
        q = (2*a2**3-9*a2*a1+27*a0)/27
        R = q**2/4+p**3/27
        Zp = poly1d([1,a2,a1,a0])
        if R <= 0: # 3 unequal real roots, VLE region
            Z = [max(Zp.roots).real,min(Zp.roots).real]
        else: # one real root for the vapor phase or liquid phase
            Z = [r.real for r in Zp.roots if r.imag == 0]
            
        return Z

class MCPR:
    """Class for multi-component Ping-Robinson equation of state.
    
    Parameters
      - s, multicomponent stream definition
      - method, default method returned by value function
      - default method for mixing rule
    
    Dependencies
      - stream dictionary from process class
      
      """
    def __init__(self, stream, method="Z", mixingRule="VDW"):
        self.method = method
        self.mixingRule = mixingRule
        self.s = stream.s
        self.mpdb = stream.mpdb
        cst = Const.Const()
        self.R = cst.R["SI"]
        self.PUnit = "Pa"
        self.TUnit = "K"
        self.vUnit = "m^3/g-mol"
        self.Pci = {}
        self.Tci = {}
        self.ai = {}
        self.bi = {}
        self.sum_yjaij = {}
        self.a = 0
        self.b = 0
        self.omegai = {}
        self.fomegai = {}
        self.s["Pressure"] = units.uc(self.s["Pressure"],\
                                      self.s["PUnit"],\
                                      self.PUnit)
        self.s["Temperature"] = units.uc(self.s["Temperature"],\
                                         self.s["TUnit"],\
                                         self.TUnit)
        self.default_kij = 0 # 0 for non-polar hydrocarbons
                             # ~ 0.5 for polar non-ideal mixtures

        for k,v in self.s["Composition"].iteritems():
            self.Pci[k] = units.uc(self.mpdb[k]["CritP"]["Value"],\
                                   self.mpdb[k]["CritP"]["VUnit"],\
                                   self.PUnit)
            self.Tci[k] = units.uc(self.mpdb[k]["CritT"]["Value"],\
                                   self.mpdb[k]["CritT"]["VUnit"],\
                                   self.TUnit)
            self.omegai[k] = self.mpdb[k]["Omega"]
            self.bi[k] = 0.07780*self.R*self.Tci[k]/(self.Pci[k])
            self.fomegai[k] = 0.37464 + 1.54226*self.omegai[k] - \
                              0.26992*self.omegai[k]**2

                    
    def value(self, *args):
        """
        """
        MCPR_Method = getattr(self, "MCPR_%s" % self.method)

        return MCPR_Method(*args)

    def mixRule(self, *args):
        """
        """
        MCPR_Method = getattr(self, "MCPR_%s" % self.mixingRule)

        return MCPR_Method(*args)

    def MCPR_VDW(self):
        """Van Der Walls mixing rule from Peng-Robinson equation of state.

        Parameters
          - 

        Returns
          -
          
          """
        
        for k,v in self.s["Composition"].iteritems():
            self.ai[k] = 0.45724*self.R**2*self.Tci[k]**2* \
                         (1+self.fomegai[k]*(1-(self.s["Temperature"]\
                                                /self.Tci[k])**0.5))**2\
                                                /(self.Pci[k])

        self.a = 0
        self.b = 0
        for ki,vi in self.s["Composition"].iteritems():
            self.b += vi*self.bi[ki]
            # Check if sum includes ii terms and ji components?
            self.sum_yjaij[ki] = 0
            for kj,vj in self.s["Composition"].iteritems():
                if ki == kj:
                    kij = 0 
                elif self.mpdb[ki].has_key("kij"):
                    kij = self.s.mpdb[ki]["kij"]
                elif self.mpdb[kj].has_key("kij"):
                    kij = self.mpdb[kj]["kij"]
                else:
                    kij = self.default_kij 
                aij = sqrt(self.ai[ki]*self.ai[kj])*(1-kij)
                self.sum_yjaij[ki] += vj*aij
                self.a += vi*vj*aij

    def MCPR_Z(self):
        """Returns the liquid and vapor phase compresibilities.

        Parameters
          -
          
        Dependencies
          - self.mixRule()
          
        Returns
          - list of gas and/or liquid compresibilities
            - for VLE: 
              [ZV, ZL]
            - for single phase liquid or gas
              [Z]
         
          """
        self.mixRule()
        P = self.s["Pressure"]
        T = self.s["Temperature"]
        B = P*self.b/self.R/T
        A = self.a*P/self.R**2/T**2
        a0 = B**3+B**2-A*B
        a1 = A-3*B**2-2*B
        a2 = B-1
        p = (3*a1-a2**2)/3
        q = (2*a2**3-9*a2*a1+27*a0)/27
        R = q**2/4+p**3/27
        Zp = poly1d([1,a2,a1,a0])
        if R <= 0: # 3 unequal real roots, VLE region
            Z = [max(Zp.roots).real,min(Zp.roots).real]
        else: # one real root for the vapor phase or liquid phase
            Z = [r.real for r in Zp.roots if r.imag == 0]
            
        return Z

    def MCPR_FugacityCoef(self, Z):
        """Fugacity coefficients from Peng-Robinson equation of state.

        Parameters
          - Z, compressibility

        Dependencies
          - self.mixRule()
          
        Returns
          - list of liquid and vapor dictionaries of fugacity
            coefficients for each coponent
          """

        self.mixRule()
        P = self.s["Pressure"]
        T = self.s["Temperature"]
        B = P*self.b/self.R/T
        sqrt2 = sqrt(2)
        fci = {}

        for ki,vi in self.s["Composition"].iteritems():
            fci[ki] = exp(self.bi[ki]/self.b*(Z-1)-log(Z-B) - \
                          (2*self.sum_yjaij[ki]/self.a-self.bi[ki]/self.b)*
                          log((Z+(1+sqrt2)*B)/ \
                              (Z+(1-sqrt2)*B))* \
                          self.a/2/sqrt2/self.b/self.R/T)
        
        return fci

if __name__ == "__main__":
    
    from mpdb.cfgMPDB import *
    from mpdb.Thermo import *
    from process import *
    import matplotlib
    matplotlib.use('TkAgg')
    from pylab import *
    import copy

    c = cfgMPDB.cfgMPDB()
    DB  = c.dbMerge(c.default_dbl)

    c1key = "7732-18-5" # water
    DB[c1key]["Omega"] = 0.344 # Currently not in any mpdb
#    c1key = "67-64-1" # acetone
    c2key = "64-17-5" # ethanol
    c3key = "67-56-1" # methanol
    T = 300.00 # K
    TUnit = "K"
    P = 101325
    PUnit = "Pa"
    
    DB[c1key]["EOS"] = pr.PCPR(DB[c1key])
    DB[c2key]["EOS"] = pr.PCPR(DB[c2key])
    DB[c3key]["EOS"] = pr.PCPR(DB[c3key])

    # default PCPR method is Pressure
    Pressure = DB[c1key]["EOS"].value(T,0.022414)
    print "%s EOS Pressure at %1.2f %s: %1.2f %s" \
          % (c1key, T, TUnit, Pressure, DB[c2key]["EOS"].PUnit)

    DB[c1key]["vp"] = vaporPres.pcvp(DB[c2key])

    def fugCoefDiff(Temp,Press,cDB):
        cDB["EOS"].method = "Z"
        Z = cDB["EOS"].value(Temp,Press)
        cDB["EOS"].method = "FugacityCoef"
        if len(Z) == 2:
            return cDB["EOS"].value(Temp,Press,Z[0]) - \
                   cDB["EOS"].value(Temp,Press,Z[1])
        else:
            return cDB["EOS"].value(Temp,Press,Z[0])

    result = optimize.fsolve(fugCoefDiff,T,args=(P, DB[c1key]))
    vpInit = units.uc(DB[c1key]['vp'].value(result),'bar', PUnit)
    
    print "Ping Robinson VP at %1.2f %s is %1.2f %s.  Antoine VP is %1.2f %s" \
          % (result,TUnit,P,PUnit, vpInit, PUnit)

    # Fit VLE data.  
    liquidComp = {c1key:0.8,
                  c2key:0.1,
                  c3key:0.1}
    liquid = process.stream(liquidComp,DB)
    liquid.s["Temperature"] = T
    liquid.s["TUnit"] = TUnit
    liquid.s["State"] = "Liquid"
    liquid.s["Basis"] = "Mole"
    liquid.s["Pressure"] = P
    liquid.s["PUnit"] = PUnit
    liquid.s["FlowRate"] = 1000
    liquid.s["FRUnit"] = "g-mole/hr"
    liquid.s["EOS"] = pr.MCPR(liquid)

    vaporComp = {c1key:0.1,
                 c2key:0.8,
                 c3key:0.1}
    vapor = process.stream(vaporComp,DB)
    vapor.s["Temperature"] = T
    vapor.s["TUnit"] = TUnit
    vapor.s["State"] = "Vapor"
    vapor.s["Basis"] = "Mole"
    vapor.s["Pressure"] = P
    vapor.s["PUnit"] = PUnit
    vapor.s["FlowRate"] = 1000
    vapor.s["FRUnit"] = "g-mole/hr"
    vapor.s["EOS"] = pr.MCPR(vapor)

    K1_meas = vapor.s["Composition"][c1key]/liquid.s["Composition"][c1key]
    print "K_%s: %1.2f" % (c1key, K1_meas)

    def dewPoint(Txin):
        vapor.s["Temperature"] = Txin[0]
        liquid.s["Temperature"] = Txin[0]
        xin = Txin[1:]

        vapor.s["EOS"].method = "Z"
        Zv = vapor.s["EOS"].value()
        vapor.s["EOS"].method = "FugacityCoef"
        vfc = vapor.s["EOS"].value(Zv[0])

        sumx = 0
        for k in mfkey:
            liquid.s["Composition"][k] = xin[mfkey.index(k)]
            sumx += xin[mfkey.index(k)]
            
        for k,v in liquid.s["Composition"].iteritems():
            if mfkey.count(k) == 0:
                liquid.s["Composition"][k] = 1-sumx
            
        liquid.s["EOS"].method = "Z"
        Zl = liquid.s["EOS"].value()
        liquid.s["EOS"].method = "FugacityCoef"
        if len(Zl) == 1:
            lfc = liquid.s["EOS"].value(Zl[0])
        else:
            lfc = liquid.s["EOS"].value(Zl[1])

        return [vapor.s["Composition"][k] - v*lfc[k]/vfc[k] for k,v in \
                liquid.s["Composition"].iteritems()]

    def bublePoint(Tyin):
        vapor.s["Temperature"] = Tyin[0]
        liquid.s["Temperature"] = Tyin[0]
        yin = Tyin[1:]

        sumy = 0
        for k in mfkey:
            vapor.s["Composition"][k] = yin[mfkey.index(k)]
            sumy += yin[mfkey.index(k)]
            
        for k,v in vapor.s["Composition"].iteritems():
            if mfkey.count(k) == 0:
                vapor.s["Composition"][k] = 1-sumy
            
        vapor.s["EOS"].method = "Z"
        Zv = vapor.s["EOS"].value()
        vapor.s["EOS"].method = "FugacityCoef"
        vfc = vapor.s["EOS"].value(Zv[0])

        liquid.s["EOS"].method = "Z"
        Zl = liquid.s["EOS"].value()
        liquid.s["EOS"].method = "FugacityCoef"
        if len(Zl) == 1:
            lfc = liquid.s["EOS"].value(Zl[0])
        else:
            lfc = liquid.s["EOS"].value(Zl[1])

        return [liquid.s["Composition"][k] - v*vfc[k]/lfc[k] for k,v in \
                vapor.s["Composition"].iteritems()]

    T_init = T
    mfkey = [c1key] 
    Tx_init = [liquid.s["Composition"][k] for k in mfkey]
    Tx_init.insert(0,T_init)
#    Txdp = optimize.fsolve(dewPoint,Tx_init)
#    print Tx_init, Txdp
    Ty_init = [vapor.s["Composition"][k] for k in mfkey]
    Ty_init.insert(0,T_init)
#    Tybp = optimize.fsolve(bublePoint,Ty_init)
#    print Ty_init, Tybp
    c.Yaws.cfg["parseOpt"] = "select"
    c.Yaws.cfg["selectFiles"] = ["VLE.csv"]
    c.Yaws.AllYawsData()
    vleDB = c.Yaws.dbl[0]
    
    vleIdx = [vleDB["VLE"].index(item) for item in vleDB["VLE"]\
              if ((c2key == item[1])) and \
              (c1key == item[5]) and \
              item[6] == 760]
    lvle = [vleDB["VLE"][idx] for idx in vleIdx \
            if vleDB["VLE"][idx][3]=='liquid']
    vvle = [vleDB["VLE"][idx] for idx in vleIdx \
            if vleDB["VLE"][idx][3]=='vapor']
    yxData = []
    for ld in lvle:
        for vd in vvle:
            if ld[7] == vd[7]:
                if ld[2] > 0 and not ld[2] == 100.0 and not ld[7] == "N/A":
                    yxData.append([vd[2],ld[2],units.uc(ld[6],"mmHg",PUnit),\
                                   ld[7]+273.15,c1key])
    
    # prune noisey data
    print yxData.pop(19), yxData.pop(18), yxData.pop(13)
    
    def resid(kij,yxData):
        result = []
        for yx in yxData:
            for k,v in liquid.s["Composition"].iteritems():
                if yx[4] == k:
                    liquid.s["Composition"][k] = yx[1]/100
                    vapor.s["Composition"][k] = yx[0]/100
                else:
                    liquid.s["Composition"][k] = 1-yx[1]/100
                    vapor.s["Composition"][k] = 1-yx[0]/100
                    
            vapor.s["Temperature"] = yx[3]
            vapor.s["Pressure"] = yx[2]
            vapor.s["EOS"].default_kij = kij
            vapor.s["EOS"].method = "Z"
            Zv = vapor.s["EOS"].value()
            vapor.s["EOS"].method = "FugacityCoef"
            vfc = vapor.s["EOS"].value(Zv[0])

            liquid.s["Temperature"] = yx[3]
            liquid.s["Pressure"] = yx[2]
            liquid.s["EOS"].method = "Z"
            liquid.s["EOS"].default_kij = kij
            Zl = liquid.s["EOS"].value()
            liquid.s["EOS"].method = "FugacityCoef"
            if len(Zl) == 1:
                lfc = liquid.s["EOS"].value(Zl[0])
            else:
                lfc = liquid.s["EOS"].value(Zl[1])
            result.append([yx[3],yx[0]/yx[1],lfc[c1key]/vfc[c1key]])
            result.append([yx[3],(1-yx[0]/100)/(1-yx[1]/100),\
                           lfc[c2key]/vfc[c2key]])
#            result.extend(yx[0]/yx[1]-lfc[c1key]/vfc[c1key])
#            result.extend((1-yx[0]/100)/(1-yx[1]/100)-lfc[c2key]/vfc[c2key])
        return result

#    kijlsq = optimize.leastsq(resid, 0.25, args=(yxData))
#    print kijlsq[0]
    kij = -0.111905585495
#    kij = 0.0182087605509
#    pd = resid(kij,yxData)
#    pT = [d[0] for d in pd]
#    pM = [d[1] for d in pd]
#    pC = [d[2] for d in pd]
#    plot(pT,pM,'ro',pT,pC,'b.')
#    xlabel('Temperature, K')
#    yLabel = 'K_%s' % c2key
#    ylabel(yLabel)
#    show()


