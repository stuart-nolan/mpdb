"""
mpdb, a basic material property database
NRTL, liquid activity coefficient model for VLE

Revision Date: 2021.08.15

SPDX-License-Identifier: BSD-2-Clause
Copyright (c) 2021 Stuart Nolan. All rights reserved.
"""
from math import exp
from scipy import zeros

class NRTL:
    def __init__(self, components={}, binaryIP={}, gasConst = 1.98588):
        """
        NRTL activity coefficient model.
    
        Parameters
           - components, dictionary of liquid components, 
             - keys are component IDs (same as in binaryIP)
             - values are fractions
           - binaryIP, dictionary of NRTL binary interaction parameters
             - binaryIP key pairs are the same keys used in components
             - ex. {('water','methanol'): (A12,A21, alpha12),
                    ('ethanol', 'water'): (A31,A13,alpha13),
                    ('ethanol', 'methanol'): (A32,A23,alpha23)} = bIP
                    {'water': 0.5, 'methanol': 0.25, 'ethanol': 0.25} = comp
             - Dechema bIPs are in cal/mol
           - gasConst, use cgs for standard Dechema bIPs
         
         Dependencies
           - scipy
        
         Returns
           - gammai as dictionary of liquid activity coefficients for components 
        
        """ 
        self.R = gasConst # kcal/mol/K for cgs system (Dechema)
        self.comp = self.comp_init(components)
        self.tau = zeros((self.nComp,self.nComp))
        self.bIP = zeros((self.nComp,self.nComp))
        self.alpha = zeros((self.nComp,self.nComp))
        self.sumj_Gjixj = zeros(self.nComp)
        self.yi = zeros(self.nComp)
        self.sumj_taujiGjixj = zeros(self.nComp)
        self.bIP_init(binaryIP)
        self.gi = zeros(self.nComp)
        
    def comp_init(self, comp):
        """
        sets up ordered and array of components keys, compenent indexs, and
        componenet mole fractions
        
        Parameters:
            comp = dictionary of coponent name keys and mole fraction values
        """
        self.nComp = len(comp)
        ci = zeros(self.nComp)
        [self.cdx,self.kdx] = zip(*enumerate(comp))
        for idx in self.cdx:
            ci[idx] = comp[self.kdx[idx]]                
        
        return ci
            
    def bIP_init(self, binaryIP):
        for idx in self.cdx:
            for jdx in self.cdx:
                if self.kdx[jdx] == self.kdx[idx]:
                    self.bIP[idx][jdx] = 0                    
                    self.alpha[jdx][idx] = 0                    
                elif (self.kdx[jdx],self.kdx[idx]) in binaryIP:
                    self.bIP[jdx][idx] = binaryIP[(self.kdx[jdx],self.kdx[idx])][0]
                    self.bIP[idx][jdx] = binaryIP[(self.kdx[jdx],self.kdx[idx])][1]                    
                    self.alpha[jdx][idx] = binaryIP[(self.kdx[jdx],self.kdx[idx])][2]
                elif (self.kdx[idx],self.kdx[jdx]) in binaryIP:
                    self.bIP[jdx][idx] = binaryIP[(self.kdx[idx],self.kdx[jdx])][1]
                    self.bIP[idx][jdx] = binaryIP[(self.kdx[idx],self.kdx[jdx])][0]                    
                    self.alpha[jdx][idx] = binaryIP[(self.kdx[idx],self.kdx[jdx])][2]
                else:
                    print("key pair: (%s, %s) not found in binaryIP") % (self.kdx[idx],self.kdx[jdx])  
            
    def gamma(self, invRT):
        comp = self.comp
        bIP = self.bIP
        alpha = self.alpha
        tau = self.tau
        sumj_Gjixj = self.sumj_Gjixj
        sumj_taujiGjixj = self.sumj_taujiGjixj
        cdx = self.cdx
        gi = self.gi
        
        for idx in cdx:
            sumj_Gjixj[idx] = 0
            sumj_taujiGjixj[idx] = 0
            for jdx in cdx:
                tau[jdx][idx] = bIP[jdx][idx] * invRT
                Gjixj = exp(-alpha[jdx][idx] * tau[jdx][idx])*comp[jdx]
                sumj_Gjixj[idx] += Gjixj
                sumj_taujiGjixj[idx] += tau[jdx][idx] * Gjixj
        
        for idx in cdx:
            sumj = 0
            for jdx in cdx:
                sumj += exp(-alpha[idx][jdx] * tau[idx][jdx]) * \
                    comp[jdx] / sumj_Gjixj[jdx] * (tau[idx][jdx] - \
                          sumj_taujiGjixj[jdx] / sumj_Gjixj[jdx])
            gi[idx] = exp(sumj_taujiGjixj[idx] / sumj_Gjixj[idx] + sumj)

        return [gi[idx] for idx in cdx]
        
def antoine_init(g, antoineConst):
    aCs = zeros((g.nComp,3))    
    for idx in g.cdx:
        aCs[idx] = antoineConst[g.kdx[idx]]

    return aCs
    
def antoine(temp,aC):
    """
    pure component pressure calculated from Antoine equation:
        log(P) = A - B/(temp + C)
    
    Parameters:
        temp - temperature in units consisten with Antoine constants 
        aC - list of Antoine equation constants [A, B, C]
    
    Returns
       - pressure in units consistent with Antoine constants 
    
    """ 
    
    return 10**(aC[0] - aC[1]/(temp+aC[2]))

def vle(temp,g,aCs,pressure):
    """
    example vle flash functions
    """
    yi = g.yi
    comp = g.comp
    invRT = 1.0/g.R/temp[0]
    gammai = g.gamma(invRT)
    
    sumy = 0
    for idx in g.cdx:
        Pi = antoine(temp[0]-273.15,aCs[idx])
        yi[idx] = comp[idx]*Pi*gammai[idx]/pressure
        sumy += yi[idx]
    
    return 1-sumy
        
from scipy.optimize import fsolve
if __name__ == "__main__":
    lComp = list(range(4))
    vComp = list(range(4))
    temp = list(range(4))
    temp[0] = 82.7
    lComp[0] = {'methanol': 0.1090, 'ethanol': 0.0810, 'water': 0}
    vComp[0] = {'methanol': 0.3550, 'ethanol': 0.2500, 'water': 0}
    temp[1] = 83.5
    lComp[1] = {'methanol': 0.0900, 'ethanol': 0.0800, 'water': 0}
    vComp[1] = {'methanol': 0.3080, 'ethanol': 0.2850, 'water': 0}
    temp[2] = 82.5
    lComp[2] = {'methanol': 0.0610, 'ethanol': 0.1420, 'water': 0}
    vComp[2] = {'methanol': 0.1950, 'ethanol': 0.3950, 'water': 0}
    temp[3] = 82.2
    lComp[3] = {'methanol': 0.1600, 'ethanol': 0.0370, 'water': 0}
    vComp[3] = {'methanol': 0.5150, 'ethanol': 0.1050, 'water': 0}
    idx = 3
    sumx = 0
    sumy = 0
    for (ki,xi) in lComp[idx].items():
        sumx += xi
        sumy += vComp[idx][ki]
        
    lComp[idx]['water'] = 1-sumx
    vComp[idx]['water'] = 1-sumy
    
    bIP = {('methanol','ethanol'): [-304.404, 460.246, 0.321], 
           ('methanol','water'): [-20.386, 766.894, 0.297],
           ('ethanol', 'water'): [-75.001, 1299.770, 0.271]}
    aC = {'methanol': [8.08097,1582.271,239.726],
          'ethanol': [8.11220,1592.864,226.184],
          'water': [8.07131,1730.630,233.426]}
    

    g = NRTL(components=lComp[idx],binaryIP=bIP)
    aCs = antoine_init(g,aC)
    vle_args = (g, aCs, 760)
    res = fsolve(vle,temp[idx]+273.15,args=vle_args)
    
    yi = {}
    yiDiff = {}
    sumy = 0
    tempDiff = temp[idx]-res[0]
    for (ki,xi) in lComp[idx].items():
        yi[ki] = g.yi[g.kdx.index(ki)]
        sumy += yi[ki]
        yiDiff[ki] = vComp[idx][ki]-yi[ki]
    
    #fit heptane toluene VLE data    
    hepTolVLE = genfromtxt("heptaneTolueneVLE.csv",delimiter=",",skip_header=1)
    Pressure = 760 #mmHg
    
    (nData, nCol) = hepTolVLE.shape
    
    bIP = {('heptane','toluene'): [-57.1265, 277.3848, 0.3010]}
    aC = {'heptane': [6.89386,1264.370,216.640],
           'toluene': [6.95007,1342.310,219.187]}
    lComp = list(range(nData))
    vComp = list(range(nData))
    temp = list(range(nData))
    gamma = list(range(nData))
    gCalc = list(range(nData))
    fugacity = list(range(nData))

    for idx in range(nData):
        lComp[idx] = {'heptane': hepTolVLE[idx,1], 'toluene': 1-hepTolVLE[idx,1]} 
        vComp[idx] = {'heptane': hepTolVLE[idx,2], 'toluene': 1-hepTolVLE[idx,2]} 
        temp[idx] = hepTolVLE[idx,0]

    gHT = NRTL(components=lComp[0],binaryIP=bIP)
    aCs = antoine_init(gHT,aC)
    cLData = [gHT.comp_init(lComp[idx]) for idx in range(nData)]
    cVData = [gHT.comp_init(vComp[idx]) for idx in range(nData)]

    for idx in range(nData):
        Pi = [antoine(temp[idx],aCs[0]),antoine(temp[idx],aCs[1])]
        gamma[idx] = [Pressure*cVData[idx][0]/(cLData[idx][0]*Pi[0]), \
                      Pressure*cVData[idx][1]/(cLData[idx][1]*Pi[1])]
    
    #gamma[0][0] = 1
    #gamma[nData-1][1] = 1

    def gammaFit(gFbIP):
        gHT.bIP_init({('heptane','toluene'): gFbIP})
        gsum = 0
        for idx in range(nData):
            gHT.comp = cLData[idx]
            invRT = 1.0/g.R/temp[idx]
            gCalc[idx] = gHT.gamma(invRT)
            gsum += ((gamma[idx][0]-gCalc[idx][0])/gamma[idx][0])**2 + \
                    ((gamma[idx][1]-gCalc[idx][1])/gamma[idx][1])**2
            
        return gsum
    
    initBIP = [-57.1265, 277.3848, 0.3010]
    from scipy.optimize import fmin
    gammaFit(initBIP)
    #finalBIP = fmin(gammaFit,initBIP,xtol=0.0000001,ftol=0.0000001)
    
