"""
mpdb: Material Property Data Base as python module
waterDielectric.py: water dielectric value,
                    Born Q, Y, X values, and
                    Debye Huckle values

Revision Date: 2021.08.15

SPDX-License-Identifier: BSD-2-Clause
Copyright (c) 2021 Stuart Nolan. All rights reserved.
"""
#from scipy.misc import derivative
#from math import isclose
import sys
class waterDielectric:
    
    def __init__(self, T=298.15, P=101325):
        """
        waterDielectric(T=298.15, P=101325)
        
        Calculates varios water dielectric and Debye Huckel values for liquid
        water (T & P range ~ 647 K > T > 250 K and 1000 MPa > P > 612 Pa)

        Dependencies:
            CoolProp - for water thermodynamic properties

        Constructor Parameters:
            T, temperature in K
            Pa, pressure in Pa 

        References:
            1 "A Formulation for the Static Permittivity of Water and Steam at 
              Temperatures from 238 K to 873 K at Pressures up to 1200 MPa, 
              Including Derivatives and Debye–Hückel Coefficients", Fernández,
              D. P.; Goodwin, A. R. H.; Lemmon, E. W.; Levelt Sengers, J.M.H.;
              Williams, R.C.;
              Journal of Physical and Chemical Reference Data 26, 1125 (1997); 
              doi: http://dx.doi.org/10.1063/1.555997
        """
        try:
            self.cP = __import__('CoolProp')
        except:
            raise ModuleNotFoundError("No module named \"CoolProp\": Try \"pip install CoolProp\"")

        self.HEOS = self.cP.AbstractState("HEOS", "Water")
        self.cPPhase = {0:'liquid',
                        1:'supercritical',
                        2:'supercritical_gas',
                        3:'supercritical_liquid',
                        4:'critical_point',
                        5:'gas',
                        6:'twophase',
                        7:'unknown',
                        8:'not_imposed'}
        self._T = T
        self.TUnit = 'K'
        self._P = P
        self.PUnit = 'Pa'

        self.rho = 0        
        self.rhoUnit = "mol/m^3"
        
        # initialize eps and Born variables
        self.eps = 0
        self.BornQ = 0
        self.BornQUnit = "Pa^-1"
        self.BornY = 0
        self.BornYUnit = "K^-1"
        self.BornX = 0
        self.BornXUnit = "K^-2"
        
        # initialize Debye Huckel values
        self.Agamma = 0
        self.Aphi= 0
        self.AV = 0
        self.AH = 0
        self.Akappa = 0
        self.AC = 0
        
        self.Tc = self.HEOS.T_critical()
        self.TcUnit = "K"
        self.Pc = self.HEOS.p_critical()
        self.PcUnit = "Pa"
        self.rhoc = self.HEOS.rhomolar_critical()
        self.rhocUnit = "mol/m^3"        
        self.fW = self.HEOS.molar_mass()
        self.fWUnit = "kg/mol"
        
        # define constants needed for calculations
        self.pi = 3.141592653589793
        # Permittivity of free space
        self.eps0 = 1.0/(4e-7*self.pi*(299792458)**2) 
        self.eps0Unit = "C^2/J/m"
        self.k = 1.380658e-23
        self.kUnit = "J/K"
        self.Na = 6.0221367e23
        self.NaUnit = "1/mol"
        self.alpha = 1.636e-40 # Mean molecular polarizability of water
        self.alphaUnit = "C^2/J/m^2"
        self.mu = 6.138e-30 # Dipole moment of water
        self.muUnit = "C*m"
        self.eC = 1.60217733e-19 # elementary charge
        self.eCUnit = "C"
        self.R = self.k*self.Na # gas constant
        self.RUnit = "J/mol/K"
        
        self.ih = [1, 1, 1, 2, 3, 3, 4, 5, 6, 7, 10, None]
        self.jh = [0.25, 1, 2.5, 1.5, 1.5, 2.5, 2, 2, 5, 0.5, 10, None]
        self.nh = [0.978224486826, -0.957771379375, 0.237511794148, 
                   0.714692244396, -0.298217036956, -0.108863472196, 
                   .949327488264e-1, -.980469816509e-2, .165167634970e-4, 
                   .937359795772e-4, -.12317921872e-9, .196096504426e-2]
        # Temperature difference used to calc d2epsdT2_P in BornX
        self.deltaT = 0.2
        self.detlaTUnit = "K"
        # Pressure difference used to calc Akappa in DebeyeHuckel
        self.deltaP = 2.000002
        self.deltaPUnit = "Pa"

    @property
    def P(self):
        return self._P
        
    @P.setter
    def P(self,P):
        self._P = P
        self.specLiqPhase()   
                        
    @property
    def T(self):
        return self._T
        
    @T.setter
    def T(self,T):
        self._T = T
        self.specLiqPhase()            

    def specLiqPhase(self):
        """
        Set phase to liquid, warn if not liquid, then call self.calc()
        """
        self.HEOS.specify_phase(self.cP.iphase_not_imposed)
        self.HEOS.update(self.cP.PT_INPUTS, self._P, self._T)
        phase = self.HEOS.phase()
        if phase == 6:
            self.HEOS.specify_phase(self.cP.iphase_liquid)
        elif phase in [0,1,3]:
            pass
        else:
            print("WARNING: P = %s %s; T = %s %s;" % (self._P,self.PUnit,self._T,self.TUnit))
            print("CoolProp phase: %s" % self.cPPhase[phase])

        self.calc()
        
    def calc(self,method='allWD'):
        """
        calc(method)
        
        method dispatcher 

        Parameters:
            method, string which may be one of:
                "eps":
                    epsEqns()
                "epsBornX":
                    epsEqns() and
                    BornXEqn()
                "epsDHEqns1":
                    epsEqns() and
                    DebyeHuckelEqns1()
                "allWD" (the default):
                    epsEqns(), 
                    BornXEqn(), 
                    DebyeHuckelEqns1(), and 
                    DebyeHuckelEqns2()

        Returns:
            calc_eqns()
        """

        calc_eqns = getattr(self, "calc_%s" % method, self.calc_allWD)

        return calc_eqns()

    def calc_allWD(self):
        """
        calc_all()
            calculates eps, Born Q, Born Y, Born X, and Debye Huckel values:
                Agamma,
                Aphi,
                AV,
                AH,
                Akappa, and
                AC
                
            Calls:
                calc_eps, 
                BornXEqn,
                DebyeHuckelEqns1, and
                DebyeHuckelEqns2
        
        i.e. the dielectric value is a function of density and T.  The property
        setters above and the calc function make it look like a function of P &
        T.
        """
        self.calc_eps()
        self.BornX = self.BornXEqn(self.eps,self.BornY,self.deltaT)

        #DebyeHuckel values require a call to epsEqns before calling calcDeyb        
        (self.Agamma, self.Aphi, 
         self.AV, self.AH) = self.DebyeHuckelEqns1(self.rho,self._T,
                                                   self.eps,self.BornQ,
                                                   self.BornY)

        (self.Akappa, self.AC) = self.DebyeHuckelEqns2(self.deltaT,
                                                       self.deltaP)
        
    def calc_eps(self):
        """calc_eps()
        
        After updateing the CoolProp HEOS instance for liquid water
        with values of P /Pa, T /K, and rho /mol/m^3, calc_eps calls
        epsEqns
        
        i.e. the dielectric value is a function of density and T.  The
        property setters above and the calc function make it look like
        a function of P & T.
        """
        # use P and T to calulate rhomolar from CoolProp HEOS instance
        self.HEOS.update(self.cP.PT_INPUTS, self._P, self._T)
        self.rho = self.HEOS.rhomolar()

        (self.eps,self.BornQ,self.BornY) = self.epsEqns(self.rho,self._T)

    def calc_epsBornX(self):
        """calc_epsBornX()
        
        After updateing the CoolProp HEOS instance for liquid water
        with values of P /Pa, T /K, and rho /mol/m^3, calc_epsBornX 
        calls:

            epsEqns and
            BornXEqn
        
        i.e. the dielectric value is a function of density and T.  The
        property setters above and the calc function make it look like
        a function of P & T.
        """

        self.calc_eps()
        
        self.BornX = self.BornXEqn(self.eps,self.BornY,self.deltaT)

    def calc_epsDHEqns1(self):
        """calc_epsDHEqns1()
        
        After updateing the CoolProp HEOS instance for liquid water
        with values of P /Pa, T /K, and rho /mol/m^3, calc_epsDHEqns1 
        calls:

            epsEqns, 
            DebyeHuckelEqns1
        
        i.e. the dielectric value is a function of density and T.  The
        property setters above and the calc function make it look like
        a function of P & T.
        """

        self.calc_eps()
        
        #DebyeHuckel values require a call to epsEqns before calling calcDeyb        
        (self.Agamma, self.Aphi, 
         self.AV, self.AH) = self.DebyeHuckelEqns1(self.rho,self._T,
                                                   self.eps,self.BornQ,
                                                   self.BornY)
    def epsEqns(self,rho,T):
        """calc_epsEqns(rho,T)
        
        water dilectric constant and Born Q and Y values as a function
        of molar density and temperature.  See also BornXEqn function
        below.
        
        Parameters:
            rho, molar density in mol/m^3
            T, temperature in K

        Returns: (eps,BornQ,BornY)
            eps, water dielectric constant at rho and T
            Born Q = d(eps)/d(P)|T/eps**2, in K^-1
            Born Y = d(eps)/d(T)|P/eps**2, in Pa^-1
            
        Notes:
            epsEqns changes self.HEOS with the arguments rho and T.
            It is the responsibility of the function calling epsEqns
            to update self.HEOS back to self.rho and self._T if that
            is what is needed.

        """
        self.HEOS.update(self.cP.DmolarT_INPUTS, rho, T)
        invTr = self.Tc/T
        rhor = rho/self.rhoc

        g = 1+self.nh[11]*rhor/(self.Tc/228.0/invTr-1)**1.2
        dgdT_rho = -1.2*self.nh[11]*rhor/(self.Tc/228.0/invTr-1)**2.2/228.0
        dgdrho_T = self.nh[11]/(self.Tc/228.0/invTr-1)**1.2/self.rhoc
        for i in range(11):
            g += self.nh[i]*rhor**self.ih[i]*invTr**self.jh[i]
            dgdT_rho += self.nh[i]*rhor**self.ih[i]*self.jh[i]* \
                        invTr**(self.jh[i]-1)*(-invTr/T)
            dgdrho_T += self.nh[i]*self.ih[i]*rhor**(self.ih[i]-1)* \
                        invTr**self.jh[i]/self.rhoc

        #calc eps
        A = self.Na*self.mu**2*rho*g/self.eps0/self.k/T
        B = self.Na*self.alpha*rho/3/self.eps0
        eps = (1+A+5*B+(9+2*A+18*B+A**2+10*A*B+9*B**2)**0.5)/4/(1-B)
        
        #for depsdP_T
        A1 = self.Na*self.mu**2*g/self.eps0/self.k/T + \
             self.Na*self.mu**2*rho/self.eps0/self.k/T*dgdrho_T
        B1 = self.Na*self.alpha/3/self.eps0
        C = 9+2*A+18*B+A**2+10*A*B+9*B**2
        depsdrho_T = 4*B1*eps/(4-4*B)+(A1+5*B1+(0.5*C**(-0.5))*(2*A1+18*B1 + \
                               2*A*A1+10*(A1*B+A*B1)+18*B*B1))/(4-4*B)
        
        #for depsdT_P
        A2 = -self.Na*self.mu**2*rho*g/self.eps0/self.k/T**2 + \
             self.Na*self.mu**2*rho/self.eps0/self.k/T*dgdT_rho
        depsdT_rho = (A2+0.5*C**(-1.0/2)*A2*(2+2*A+10*B))/(4-4*B)        

        #Born Q and Y functions
        drhodP_T = self.HEOS.first_partial_deriv(self.cP.iDmolar,self.cP.iP,
                                                 self.cP.iT)
        dPdT_rho = self.HEOS.first_partial_deriv(self.cP.iP,self.cP.iT,
                                                 self.cP.iDmolar)
        invEps2 = 1.0/eps**2
        Q=invEps2*depsdrho_T*drhodP_T #depsdP_T/eps**2
        Y=invEps2*(depsdT_rho - depsdrho_T*dPdT_rho*drhodP_T) #depsdT_P/eps**2
        
        return (eps,Q,Y)

    def BornXEqn(self,eps,Y,deltaT=0.1):
        """
        calcBornX(eps,Y,deltaT=0.1)
        
        calculate Born X = d2(eps)/d(T)2|P/eps**2 - 2.0*eps*Y**2, in K^-2
        
        Parameters:
            eps, molar density in mol/m^3
            Y, temperature in K
            deltaT, temperature difference in K for d2(eps)/d(T)2|P derivative 

        Returns: 
            BornX
            
        Notes:
            epsEqns must be called before calling calcBornX                       
        """
        TH = self._T+deltaT
        TL = self._T-deltaT

        self.HEOS.update(self.cP.PT_INPUTS, self._P, TH)
        (epsH,QH,YH) = self.epsEqns(self.HEOS.rhomolar(),TH)
        
        self.HEOS.update(self.cP.PT_INPUTS, self._P, TL)
        (epsL,QL,YL) = self.epsEqns(self.HEOS.rhomolar(),TL)
        d2epsdT2_P = (YH*epsH**2-YL*epsL**2)/2/deltaT
        
        self.HEOS.update(self.cP.DmolarT_INPUTS, self.rho, self._T)
        return d2epsdT2_P/eps**2 - 2.0*eps*Y**2
    
    def DebyeHuckelEqns1(self,rho,T,eps,Q,Y):
        """
        DebyeHuckelEqns1(rho,T,eps,Q,Y)
        
        Debye Huckle coefficient: Agamma = (2*Pi*Na*rho*fW)**(1/2)*[e**2/
                                           (4*Pi*eps*eps0*k*T)]**3/2 
                                           in (kg/mol)^(1/2)
        osmotic coefficient: Aphi = Agamma/3 in (kg/mol)^(1/2)
        apparent molar volume: AV = -4*R*T*(dAphi/dP)|T
                                    in cm^3*kg^(1/2)*mol^(-3/2)
        apparent molar enthalpy: AH = 4*R*T*(dAphi/dT)|P 
                                      in (kg/mol)^(1/2)*Pa*m^3/mol
        
        Parameters:
            rho, molar density in mol/m^3
            T, temperature in K
            eps, water dielectric constant
            Born Q = d(eps)/d(T)|P/eps**2, in K^-1
            Born Y = d(eps)/d(P)|T/eps**2, in Pa^-1

        Returns: 
            (Agamma, Aphi, AV, AH)
        """

        Agamma = (2*self.pi*self.Na*rho*self.fW)**(1.0/2.0)* \
                 (self.eC**2/(4*self.pi*eps*self.eps0*self.k*T))**(3.0/2.0)

        Aphi = Agamma/3

        drhodP_T = self.HEOS.first_partial_deriv(self.cP.iDmolar,self.cP.iP,
                                                 self.cP.iT)
        AV = 2*Aphi*self.R*T*(3*Q*eps-drhodP_T/rho)

        drhodT_P = self.HEOS.first_partial_deriv(self.cP.iDmolar,self.cP.iT,
                                                 self.cP.iP)
        AH = -6*Aphi*self.R*T*(1+T*Y*eps-T*drhodT_P/3/rho)
        
        return (Agamma, Aphi, AV, AH)

    def DHEqns1_TP(self,P,T):
        self.HEOS.update(self.cP.PT_INPUTS, P, T)
        rhoD = self.HEOS.rhomolar()
        (eps,Q,Y) = self.epsEqns(rhoD,T)
        return self.DebyeHuckelEqns1(rhoD,T,eps,Q,Y)
        
    def DebyeHuckelEqns2(self,deltaT,deltaP):
        """
        DebyeHuckelEqns2(deltaT,deltaP)
        
        apparent molar compressibility: 
          Akappa = (dAV/dP)|T in m^3*kg^(1/2)*mol^(-3/2)*Pa^-1
        apparent molar heat capacity: 
          AC = (dH/dT)|P in (kg/mol)^(1/2)*Pa*m^3/mol/K
        
        Parameters:
            deltaT, temperature difference in K for d2(eps)/d(T)2|P derivative 
            deltaP, pressure difference in Pa for d2(eps)/d(T)2|P derivative

        Returns: 
            (Akappa, AC)
        """
        
        #self.HEOS.update(self.cP.PT_INPUTS, self._P+deltaP, self._T)
        #rhoD = self.HEOS.rhomolar()
        #(epsPH,QPH,YPH) = self.epsEqns(rhoD,self._T)
        (AgammaPH, AphiPH, AVPH, AHPH) = self.DHEqns1_TP(self._P+deltaP,
                                                         self._T)
        (Agamma2PH, Aphi2PH, AV2PH, AH2PH) = self.DHEqns1_TP(self._P+2*deltaP,
                                                             self._T)
        (AgammaPL, AphiPL, AVPL, AHPL) = self.DHEqns1_TP(self._P-deltaP,
                                                         self._T)
        (Agamma2PL, Aphi2PL, AV2PL, AH2PL) = self.DHEqns1_TP(self._P-2*deltaP,
                                                             self._T)
        Akappa = (AV2PL - 8*AVPL + 8*AVPH - AV2PH)/12/deltaP

        (AgammaTH, AphiTH, AVTH, AHTH) = self.DHEqns1_TP(self._P,
                                                         self._T+deltaT)
        
        (Agamma2TH, Aphi2TH, AV2TH, AH2TH) = self.DHEqns1_TP(self._P,
                                                             self._T+2*deltaT)
        
        (AgammaTL, AphiTL, AVTL, AHTL) = self.DHEqns1_TP(self._P,
                                                         self._T-deltaT)
        (Agamma2TL, Aphi2TL, AV2TL, AH2TL) = self.DHEqns1_TP(self._P,
                                                             self._T-2*deltaT)
        AC = (AH2TL - 8*AHTL + 8*AHTH - AH2TH)/12/deltaT
        
        self.HEOS.update(self.cP.DmolarT_INPUTS, self.rho, self._T)
        return (Akappa,AC)
    
def validate(wDC,HAVE_TABULATE):
    """
    Selected results from Ref 1. for comparison.
    """
        
    header_r1T17 = ["T /K", "P /Pa", "Aphi\n/(kg/mol)^(1/2)",
                    "AV\n/m^3*kg^(1/2)/mol^(3/2)", "AH/R/T\n/(kg/mol)^(1/2)",
                    "Akappa\n/m^3*kg^(1/2)/mol^(3/2)/Pa",
                    "AC/R\n/(kg/mol)^(1/2)"]
    # Reference 1, Table 17
    r1T17 = [[300,101325,0.39251,1.9275e-6,0.81448,-4.3409e-15,3.9014],
             [300,100000000,0.37505,1.5854e-6,0.72666,-2.72e-15,3.1524],
             [500,10000000,0.67109,2.0023e-5,4.9811,-2.1966e-13,29.970],
             [550,1000000000,0.41568,2.4520e-6,1.3084,-2.1672e-15,2.8536],
             [600,100000000,0.80562,3.2686e-5,6.5128,-3.3565e-13,36.568],
             [800,100000000,2.0669,6.1662e-4,38.835,-2.0833e-10,239.33],
             [800,1000000000,0.56392,6.4374e-6,1.8592,-8.3153e-15,3.4695]]

    valData_r1T17 = []
    ff = lambda n : float("{:.5}".format(n))
    for row in r1T17:
        wDC._T = row[0]
        wDC.P = row[1]
        valData_r1T17.append([wDC.T, wDC.P, ff(wDC.Aphi), ff(wDC.AV),
                              ff(wDC.AH/wDC.R/wDC.T),ff(wDC.Akappa),
                              ff(wDC.AC/wDC.R)])
    if HAVE_TABULATE:
        print("Selected values from Fernández et. al 1997 (Ref. 1), Table 17:")
        print(tabulate(r1T17,headers=header_r1T17))
        print("\n Results from %s:" % __file__)
        print(tabulate(valData_r1T17,headers=header_r1T17))

    return (header_r1T17,r1T17,valData_r1T17)

if __name__ == "__main__":
    try:
        from tabulate import tabulate
        HAVE_TABULATE = True
    except:
        HAVE_TABULATE = False

    wDC = waterDielectric(T=300,P=101325)
    (h_r1T17,rD_r1T17,vD_r1T17) = validate(wDC,HAVE_TABULATE)
