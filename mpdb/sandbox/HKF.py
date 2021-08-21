"""
mpdb: Material Property Data Base as python module
HKF.py: aqueous electrolyte EOS 

Revision Date: 2021.08.15

SPDX-License-Identifier: BSD-2-Clause
Copyright (c) 2021 Stuart Nolan. All rights reserved.
"""

from mpdb.waterDielectric import waterDielectric as wD
from math import log

class aqEOS(wD):
    def __init__(self,T=298.15,P=101325):
        """Calculate standard partial molal volume, heat capacity, Gibbs
        free energy, and entropy of aqueous metal complexes (e.g. ions,
        associated salts, etc.) between 0-1000 deg. C and 1-5000
        bar. Standard partial molal quantites are based on molality
        (mole component per kg of water) concentration scale.
        
        Reference 1. below refers to the equations used in this module
        as the "revised HKF equation of state."  Reference 2 and the
        code for supcrt92 is used to implement selected equations in
        reference 1.  As Ref. 1. does not give an equation for
        standard partial molal entropy

        Dependencies: - mpdb.waterDielectric, which in turn depends
        on: - CoolProp

        References:

        1. "Prediction of the thermodynamic properties of aqueous metal 
            complexes to 1000°C and 5 kb"
           D. A. Sverjensky, E. L. Shock, and H. C. Helgeson,
           Geochimica et Cosmochimica Acta, Vol. 61, No. 7, pp. 1359-1412, 
           1997
            
        2. @article{JOHNSON1992899,
           title = "SUPCRT92: A software package for calculating the standard 
           molal thermodynamic properties of minerals, gases, aqueous 
           and reactions from 1 to 5000 bar and 0 to 1000°C",
           species journal = "Computers & Geosciences",
           volume = "18",
           number = "7",
           pages = "899 - 947",
           year = "1992",
           note = "",
           issn = "0098-3004",
           doi = "http://dx.doi.org/10.1016/0098-3004(92)90029-Q",
           url = "http://www.sciencedirect.com/science/article/pii/009830049290029Q",
           author = "James W. Johnson and Eric H. Oelkers and Harold C. Helgeson",
           keywords = "SUPCRT92",
           keywords = "Equations of state",
           keywords = "Standard molal thermodynamic properties",
           keywords = "Chemical equilibrium",
           keywords = "Minerals",
           keywords = "Gases",
           keywords = "Aqueous species",
           keywords = "HO",
           keywords = "Thermodynamics",
           keywords = "Geochemistry"
           }            

        Notes: 
          In supcrt92:
          - Z in aqsps is -1/eps (i.e. Born92)
          - Z in omeg92 is ion charge
          - gref = 0 when reference conditions are 298.15 K and 1 bar
          - The magnitude for EOS coefficients and omega is adjusted:
            a(1,i)  = a(1,i)*1.0d-1
            a(2,i)  = a(2,i)*1.0d2
            a(4,i)  = a(4,i)*1.0d4
            c(2,i)  = c(2,i)*1.0d4
            wref(i) = wref(i)*1.0d5

        """
        wD.__init__(self,T=T,P=P)

        self._Tref = 298.15 # K
        self._Pref = 100000 # Pa
        self.gref = 0
        self.YPrTr = 0
        self.ZPrTr = 0
        self.calc('waterRefDC')
        
        self.c = [-0.2037662E+01, 0.5747000E-02, -0.6557892E-05,
                  0.6107361E+01, -0.1074377E-01, 0.1268348E-04]
        self.cc = [0.3666666E+02, -0.1504956E-9, 0.5017997E-13]

        self.g = 0
        self.dgdP = 0
        self.dgdT = 0
        self.d2gdT2 = 0

        self.eta = 1.66027e5*4.1868 #Angstroms*J/mol
        self.psi = 260000000 # Pa
        self.theta = 228 # K
    
    # overide P & T properties from waterDielectric class
    @property
    def P(self):
        return self._P
        
    @P.setter
    def P(self,P):
        self._P = P
        self.calc('gFun')
                        
    @property
    def T(self):
        return self._T
        
    @T.setter
    def T(self,T):
        self._T = T
        self.calc('gFun')

    @property
    def Pref(self):
        return self._Pref
        
    @Pref.setter
    def Pref(self,P):
        self._Pref = P
        self.calc('waterRefDC')
                        
    @property
    def Tref(self):
        return self._Tref
        
    @Tref.setter
    def Tref(self,T):
        self._Tref = T
        self.calc('waterRefDC')

    def calc_waterRefDC(self):
        """
        calc_waterRefDC()
        
        Call once upon initialization.  Currently only supports:
            Tref = 298.15 and 
            Pref = 100000 Pa

        Computes:
            YPrTr, Born Y function at Pref and Tref
            ZPrTr, -1/eps at Pref and Tref

        Notes:
              from supcrt92 in SUP92D.F lines 1538-1545:
                  assignment of [Z,Y]PrTr to Johnson-Norton (1991) 
                  values for distribution version

                  ZPrTr = -0.1278034682d-1
                  YPrTr = -0.5798650444d-4 

              where:
                  Tref = 298.15 K
                  Pref = 1 bar
        """
        P = self._P
        T = self._T
        self._P = self._Pref
        self._T = self._Tref
        self.calc_eps()
        self.YPrTr = self.BornY
        self.ZPrTr = -1.0/self.eps
        #self.gref = ?? - how to compute/implement new ref T & P?
        
        #reset eps & HEOS to T & P        
        self._P = P
        self._T = T
        self.calc_eps()

    def calc_gFun(self):
        """
        calc_gFun()
        
        Computes self.g, self.dgdP, self.dgdT, and self.d2gdT2 at T & P 
        using equations given by Shock et al. (1991) 
            
        Parameters:
                
        Dependencies:
            
        Computes: 
            self.g, a solvent function of temperature and density in Angstroms 
                    given by Shock et al. (1992)
            self.dgdP, in Angstrom/Pa
            self.dgdT, in Angstrom/K
            self.d2gdT2, in Angstrom/K**2
            
        NOTES:
            gFun is a translation of gShok2 from REAC92D.F that is a
            part of supcrt92.  In gShok2, the units are reported as:
                TdegC.............. deg C
                Pbar............... bar
                D ................. g/cm**3
                beta, dgdP ........ bar**(-1)
                alpha, dgdT ....... K**(-1)
                daldT, d2gdT2 ..... K**(-2)

            The supcrt92 g likely has units of Angstroms as supcrt92 does not 
            redimensionalize g before it is added to dimensional values such as
            reref.

            alpha = beta*(dP/dT)|rho = coefficient of isobaric thermal expansion
            daldT = beta*(d2P/dT2)|rho
            alpha = HEOS.isobaric_expansion_coefficient() in K**-1
            beta = 1/D*(dD/dP)|T = coefficient of isothermal compressibility
            beta = HEOS.isothermal_compressibility() in Pa**-1
        """

        self.calc_epsBornX()

        self.Z = -1.0/self.eps
        alpha = self.HEOS.isobaric_expansion_coefficient() # K**-1
        beta = self.HEOS.isothermal_compressibility() # Pa**-1
        daldT = beta*self.HEOS.second_partial_deriv(self.cP.iP, self.cP.iT, 
                                                    self.cP.iDmass, self.cP.iT,
                                                    self.cP.iDmass) # K**-2
        beta = beta*100000.0 # bar**-1
        
        D = self.HEOS.rhomass()/1000.0 # g/cc

        TdegC = self._T-273.15
        Pbar = self._P/100000.0

        if D >= 1.0 or TdegC > 1000.0 or Pbar > 5000.0: return
                
        a = self.c[0] + self.c[1]*TdegC + self.c[2]*TdegC**2
        b = self.c[3] + self.c[4]*TdegC + self.c[5]*TdegC**2

        OnemD = 1.0 - D
        lnOnemD = log(OnemD)
        self.g = a*OnemD**b

        OnemDpbm1 = OnemD**(b - 1.0)
        OnemDpbm2 = OnemD**(b - 2.0)
        dgdD   = -a * b * OnemDpbm1

        dadT   =   self.c[1] + 2.0*self.c[2]*TdegC
        dadTT  =   2.0 * self.c[3]
        dbdT   =   self.c[4] + 2.0*self.c[5]*TdegC
        dbdTT  =   2.0 * self.c[5] 
        
        dDdT   = - D * alpha
        dDdP   =   D * beta
        dDdTT  = - D * (daldT - alpha**2)

        Db     = OnemD**b

        dDbdT  = -b * OnemDpbm1 * dDdT + lnOnemD * Db  * dbdT

        dDbdTT = -1*(b * OnemDpbm1 * dDdTT + \
                 OnemDpbm1 * dDdT * dbdT + b * dDdT * \
                 (-(b-1.0) * OnemDpbm2 * dDdT + \
                 lnOnemD * OnemDpbm1 * dbdT)) + \
                 lnOnemD * OnemD**b * dbdTT - \
                 OnemD**b * dbdT * dDdT / OnemD + \
                 lnOnemD * dbdT * dDbdT

        self.dgdP   = dgdD * dDdP
        self.dgdT   = a*dDbdT + Db*dadT
        self.d2gdT2 = a*dDbdTT + 2.0*dDbdT*dadT + Db*dadTT
        
        if ((TdegC < 155.0) or (Pbar > 1000.0) or (TdegC > 355.0)): return
        
        kmPbar = 1000.0 - Pbar
        kmPbar2 = kmPbar*kmPbar
        TdegCm155d300 = (TdegC - 155.0)/300.0
        
        ft     = (TdegCm155d300)**4.8 + \
                 self.cc[0]*(TdegCm155d300)**16
        
        dftdT  = 4.8/300.0*(TdegCm155d300)**3.8 + \
                 16.0/300.0*self.cc[0]*(TdegCm155d300)**15
        
        dftdTT = 3.8*4.8/300.0**2*(TdegCm155d300)**2.8 + \
                 15.0*16.0/300.0**2*self.cc[0]*(TdegCm155d300)**14
        
        fp     = self.cc[1]*kmPbar2*kmPbar + self.cc[2]*kmPbar2*kmPbar2
        
        dfpdP  = -3.0*self.cc[1]*kmPbar2 - \
                   4.0*self.cc[2]*kmPbar2*kmPbar
        
        f      = ft * fp
        dfdP   = ft * dfpdP
        dfdT   = fp * dftdT
        d2fdT2 = fp * dftdTT
        
        self.g      = self.g      - f
        self.dgdP   = self.dgdP   - dfdP
        self.dgdT   = self.dgdT   - dfdT
        self.d2gdT2 = self.d2gdT2 - d2fdT2

    def omega(self,charge,name,wref):
        """
        omega(charge,name,wref)
        
        Computes the conventinal Born coefficient omega (w) of the 
        current aqueous species, dwdP, dwdP, and dw2dT2 as a 
        function of g, dgdP, dgdT, d2gdT2, wref, and Z using
        equations given by Johnson et al. (1991).

        Parameters:
            charge, interger charge of component name
            name, string ID of component
            wref, ref omega found in slopYY.dat for each component
            
        Denpendencies:
            call calc_gFun() at least once or for every T &/or P change 
            before calling omega
            
        Returns: (w,dwdP,dwdT,d2wdT2)
            w, conventional Born coefficient in J/mol
            dwdP, d(w)/d(P)|T in J/mol/Pa
            dwdT, d(w)/d(T)|P in J/mol/K
            d2wdT2, d2(w)/d(T)2|P in J/mol/K**2

        Notes:
            3.082 is the effective electrostatic radius of the 
            hydrogen ion at 1 bar and 298.15 K in Angstrom
        """
        
        if charge == 0 or name == 'H+': #neutral aqueous species or H+ 
            omega  = wref
            dwdP   = 0.0
            dwdT   = 0.0
            d2wdT2 = 0.0

            return (omega,dwdP,dwdT,d2wdT2)

        else: #charged aqueous species other than H+
            chrg2 = charge*charge
            Hionpg = 3.082 + self.g
            Hionpg2 = Hionpg*Hionpg
            
            reref = chrg2 / (wref/self.eta + charge/(3.082+self.gref))
            re = reref + abs(charge) * self.g
            re2 = re*re
            omega  = self.eta * (chrg2/re - charge/Hionpg)
            Z3 = abs(chrg2*charge)/re2 - charge/Hionpg2
            Z4 = chrg2*chrg2/(re2*re) - charge/(Hionpg2*Hionpg)
            dwdP   = -self.eta * Z3 * self.dgdP
            dwdT   = -self.eta * Z3 * self.dgdT
            d2wdT2 = 2*self.eta*Z4*self.dgdT*self.dgdT - self.eta*Z3*self.d2gdT2

            return (omega,dwdP,dwdT,d2wdT2)

    def aqHKF(self,comp):
        """
        aqHKF(comp)
        
        Computes:

        Parameters:
            comp, 1 component from mpdb compatible slop dictionary

        Denpendencies:
            calc_gFun() (only if T &/or P change) for omega(Z,name,wref)
            
        Returns: (Gaqs,Haqs,Cpaqs,Saqs,Vaqs)
            Gaqs: Gibbs free energy of formation in J/mol
            Haqs: enthalpy of formation in J/mol
            Cpaqs: heat capacity in J/mol/K
            Saqs: entropy of formation in J/mol/K
            Vaqs: volume in m^3/mol

        Notes:

        """
        (omega,dwdP,dwdT,d2wdT2) = self.omega(comp['charge'], comp['name'],
                                              comp['omega'])

        TmTheta = self._T-self.theta
        invTmTheta = 1.0/TmTheta
        invTmTheta2 = invTmTheta*invTmTheta
        TmTref = self._T-self._Tref
        invTrefmTheta = 1.0/(self._Tref-self.theta)
        PmPref = self._P-self._Pref
        psipP = (self.psi+self._P)
        psipPref = self.psi+self._Pref
        lnpsipPdpsipPref = log(psipP/psipPref) 
        mZm1 = -self.Z - 1.0
        mZPrTrm1 = -self.ZPrTr - 1.0
        TdTref = self._T/self._Tref
        
        #Ref 1) and 2) Vaqs eqn
        #note 1 J/(mol*Pa) = 1 m**3/mol
        VQterm = (-omega*self.BornQ+mZm1*dwdP)
        Vaqs = comp['a1'] + comp['a2']/psipP + invTmTheta*(comp['a3'] + \
               comp['a4']/psipP) + VQterm

        #Ref 2) Saqs eqn - no ref 1) eqn for Saqs
        SYterm = omega*self.BornY - mZm1*dwdT - comp['omega']*self.YPrTr
        Saqs = comp['eFS298K'] + comp['c1']*log(TdTref) - \
               comp['c2']/self.theta*(invTmTheta - \
                   invTrefmTheta + (1.0/self.theta) * \
                   log(invTrefmTheta*TmTheta/TdTref)) + \
               (comp['a3']*PmPref + \
                comp['a4']*lnpsipPdpsipPref) * invTmTheta2 + SYterm
        
        #Ref 1) & 2) Cpaqs:
        CpXtrm = omega*self._T*self.BornX + 2.0*self._T*self.BornY*dwdT - \
                 self._T*mZm1*d2wdT2
        Cpaqs  = comp['c1'] + comp['c2']*invTmTheta2 - \
                 (2.0*self._T*invTmTheta2*invTmTheta) * \
                 (comp['a3']*PmPref + comp['a4']*lnpsipPdpsipPref) + CpXtrm

        #Ref 2) Haqs - no ref 1) eqn for Haqs
        HYterm = omega*mZm1 + omega*self._T*self.BornY - \
                 self._T*mZm1*dwdT - \
                 comp['omega']*mZPrTrm1 - \
                 comp['omega']*self._Tref*self.YPrTr
        Haqs = comp['eFH298K']*1000 + comp['c1']*TmTref - \
               comp['c2']*(invTmTheta - invTrefmTheta) + \
               comp['a1']*PmPref + comp['a2']*lnpsipPdpsipPref + \
               (comp['a3']*PmPref + comp['a4']*lnpsipPdpsipPref) * \
               ((2.0*self._T - self.theta)*invTmTheta2) + \
               HYterm    
        
        #Ref 2) Gaqs eqn
        GZterm = omega*mZm1 - comp['omega']*mZPrTrm1 + \
                 comp['omega']*self.YPrTr*TmTref
        Gaqs = comp['eFG298K']*1000 - comp['eFS298K']*TmTref - \
               comp['c1']*(self._T*log(TdTref)-self._T+self._Tref) + \
               comp['a1']*PmPref + comp['a2']*lnpsipPdpsipPref - \
               comp['c2']*((invTmTheta - invTrefmTheta) * \
                   ((self.theta-self._T)/self.theta)-self._T/self.theta**2 * \
                   log((self._Tref*TmTheta)/self._T*invTrefmTheta)) + \
               invTmTheta * (comp['a3']*PmPref + \
               comp['a4']*lnpsipPdpsipPref) + GZterm
        """
        #Ref 1) Gaqs eqn
        GZterm = comp['omega']*(-self.Z + self.ZPrTr + self.YPrTr*TmTref) + \
                 (omega - comp['omega'])*mZm1
        
        Gaqs = comp['eFG298K']*1000 - comp['eFS298K']*TmTref - \
               comp['c1']*(self._T*log(TdTref)-TmTref) + \
               comp['a1']*PmPref + comp['a2']*lnpsipPdpsipPref - \
               comp['c2']*((invTmTheta - invTrefmTheta) * \
                   (1.0-self._T/self.theta)-self._T/self.theta**2 * \
                   log(invTrefmTheta*TmTheta/TdTref)) + \
               comp['a3']*log(PmPref*invTmTheta) + \
               comp['a4']*lnpsipPdpsipPref*invTmTheta + GZterm
        """               
        return (Gaqs,Haqs,Cpaqs,Saqs,Vaqs)

def llnlK(T,A):
    """
    return log(K,10) as a function of T given parameter list A.
    
    """
    #T = T-273.15
    return A[0]+A[1]*T+A[2]/T+A[3]*log(T,10)+A[4]/T/T +A[5]*T*T

if __name__ == "__main__":
    from mpdb.utils import dbLoad
    import os
    from math import exp
    location = os.path.dirname(os.path.realpath(os.path.abspath(__file__)))
    dataDir = os.path.join(location,'data')
    (mpdb, mpdbmd) = dbLoad(interactive=True,dataDir=dataDir)
    tEOS = aqEOS()
    tEOS.P = 101325 #Pa
    tEOS.T = 300+273.15 #K
    # get pressue of sat. liquid water (Q = 0) at T
    # lln log10K's from supcrt92 are along saturation line
    tEOS.HEOS.update(tEOS.cP.QT_INPUTS,0,tEOS.T)
    # but waterDielectric implementation does not yet support saturation
    # conditions
    tEOS.P = tEOS.HEOS.p()+9
    """
    aqComp = [('HS-',-1),('H+',-1),('H2S,AQ',1)]
    aqComp = [('H+',-1),('CO3-2',-1),('HCO3-',1)]
    A_NH4 = [-1.4527e+001,-5.0518e-003,3.0447e+003,6.0865,4.7515e+001, 0]
    A_NaCl = [159.605, 0.084294, -3975.6, -66.857, 0, -4.9364e-05]
    A_HCO3 = [107.8975, 0.03252849, -5151.79, -38.92561, 563713.9, 0]
    A_H2S = [-782.43945, -0.361261, 20565.7315, 328.67496, 0, 1.6722e-4]
    """
    aqComp = [('NH3,AQ',-1),('H+',-1),('NH4+',1)]
    llnlKP = [-1.4527e+001,-5.0518e-003,3.0447e+003,6.0865,4.7515e+001, 0] 
    sComp = []
    dGVals = []
    for (name,pref) in aqComp:
        (Gaqs,Haqs,Cpaqs,Saqs,Vaqs) = tEOS.aqHKF(mpdb[name])
        #Gaqs = mpdb[name]['eFG298K']*1000 #J/mol
        dGVals.append((pref,Gaqs))
    #dGVals.append((sComp[0][0],sComp[0][1],mpdb[sComp[0][0]]['eFG298K']))
    dGr = 0
    for (prefix,dGc) in dGVals:
        dGr = dGr + prefix*dGc
    lnKr = -dGr/tEOS.R/tEOS.T #log base e
    """
    iphreeqc databases report log_k as log10(K) = -deltaG/R/298.15
    """
    #print lnKr in log base 10
    print("HKF log10(K) @ T = %s: %s" % (tEOS.T,lnKr/log(10))) 
    print("llnl log10(K) @ T = %s: %s" % (tEOS.T,llnlK(tEOS.T,llnlKP))) 
    #CuCl+ at Psat & T = 250 deg C lnKr/2.3025851 = -2.19 (Table 13, pg 1401 ref 1)
