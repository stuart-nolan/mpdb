"""
SPDX-License-Identifier: BSD-2-Clause
Copyright (c) 2021 Stuart Nolan. All rights reserved.
"""
from math import sqrt, exp
import os
import re
import csv
import inspect
    
def fullerADV(formula, ring=False, numRing = 1):
    """
    fullerADV(formula, ring=False, numRing = 1)
    
    Atomic Diffusion Volume for fullerD() function.  Implements calculations
    from Table 11-1 pg 11.11 of Ref. below.
    
    Refs: 
        "THE PROPERTIES OF GASES AND LIQUIDS", Bruce E. Poling, John M. 
        Prausnitz, John P. O’Connell, 5th ed., McGRAW-HILL, NY, (2001), 
        pg 11.11.
        
        Fuller, E. N., K. Ensley, and J. C. Giddings: J. Phys. Chem., 73: 3679 
        (1969).

    Parameters:
        formula, str elements that make up a component e.g. H2O, CH4, CO2, etc.
        ring, boolean, True if the compond contains an aromatic or heterocyclic
              ring
        numRing, int, number of heterocyclic or aromatic rings in the structure
        
        NOTES: 
            if an element in the formula is not found in the dictionary 
            atomic [C, H, O, N, Cl, S, I, Br, F], this functions skips that 
            element with NO notification...
            
            it is unclear to the author of this function what to do about a 
            compond with more than one aromatic or heterocyclic ring...
            
            numRing can be used to multiply the atomic diffusion volume for a
            ring by the number of rings in the compond.
        
    Returns:
        Atomic diffusion volume for the fullerD() function
    """
    defined = {'H2':6.12,'D2':6.84,'He':2.67,'N2':18.5,'O2':16.3,'Ne':5.98,
               'air':19.7,'Ar':16.2,'Kr':24.5,'Xe':32.7,'CO':18.0,'CO2':26.9,
               'N2O':35.9,'H3N':20.7,'H2O':13.1,'F6S':71.3,
               'CI2':38.4,'Br2':69.0,'SO2':41.8}

    if formula in defined.keys():
        return defined[formula]

    reElem = re.compile(r'([A-Z][a-z]*)(\d*)')
    aDV = 0
    atomic = {'C':15.9,'H':2.31,'O':6.11,'N':4.54,'Cl':21.0,'S':22.9,'I':29.8,
              'Br':21.9,'F':14.7}
              
    for (elem,count) in reElem.findall(formula):
        try:
            count = float(count)
        except:
            count = 1
        if elem in atomic.keys():
            aDV += atomic[elem]*count
    if ring:
        aDV += -18.3*numRing
        
    return aDV

def fullerD(T, P, mWA, aDVA, mWB, aDVB):
    """
    fullerD(T, P, mWA, aDVA, mWB, aDVB)
    
    Binary vapor diffusion coefficient correlation at low pressure (ideal gas
    or P < 30 bar) a la Fuller et al.  
    
    Refs:
        "THE PROPERTIES OF GASES AND LIQUIDS", Bruce E. Poling, John M. 
        Prausnitz, John P. O’Connell, 5th ed., McGRAW-HILL, NY, (2001), 
        pgs 11.10-12.
        
        Fuller, E. N., and J. C. Giddings: J. Gas Chromatogr., 3: 222 (1965).
        
        Fuller, E. N., P. D. Schettler, and J. C. Giddings: Ind. Eng. Chem., 
        58(5): 18 (1966).
        
        Fuller, E. N., K. Ensley, and J. C. Giddings: J. Phys. Chem., 73: 3679 
        (1969).
    
    Parameters:
        T, temperature in K
        P, pressure in Pa
        mWA, molecular weight of component A
        aDVA, atomic diffusion volume for component A
        mWB, molecular weight of component B
        aDVB, atomic diffusion volume for component B
    
        fullerADV(formula,ring,numRing) calculates the atomic diffusion volume
    
    Returns:
        Diffusion coefficient of A in B in cm^/s
    """
    return 101.3*T**(1.75)*(1/mWA+1/mWB)**(0.5)/P/(aDVA**(1/3)+aDVB**(1/3))**2
    
def birdD(T, P, mWA, PcA, TcA, mWB, PcB, TcB, a=2.745E-4,b=1.823):
    """
    Binary vapor diffusion coefficient correlation at low pressure (ideal gas
    or P < 30 bar)
    
    for nonpolar vapor pairs (e.g. component A or B is not H2O) or
    for water and a nonpolar vapor
    
    Refs:
        BSL pg 505, eqn 16.3-1
        J. C. Slattery and R. B. Bird, A.I.Ch.E. Journal, 4, 137-142 (1958)

    Parameters:
        T, temperature in K
        P, pressure in atm
        mWA, molecular weight of component A
        PcA, critical pressure of component A in atm
        TcA, critical temperature of component A in K
        mWB, molecular weight of component B
        PcB, critical pressure of component B in atm
        TcB, critical temperature of component B in K
        a,b, a=3.640E-4 and b=2.334 if component A or B is water,
             a=2.745E-4 and b=1.823 if neither component A or B is water
        
    Returns:
        Diffusion coefficient of A in B in cm^/s
    
    """
    
    return (PcA*PcB)**(1/3)/P*(TcA*TcB)**(5/12)*(1/mWA+1/mWB)**(1/2)*a*(T/sqrt(TcA*TcB))**(b)

def chapEnskD(T, P, mWA, sigmaA, epsilonA, mWB, sigmaB, epsilonB):
    """
    Binary dilute vapor diffusion coefficients at low pressure (ideal gas or 
    P < 30 bar) for nonpolar, spherical, monatomic molecules a la 
    Chapman Enskog et al.  
    
    Refs:
        "THE PROPERTIES OF GASES AND LIQUIDS", Bruce E. Poling, John M. 
        Prausnitz, John P. O’Connell, 5th ed., McGRAW-HILL, NY, (2001), 
        pgs 11.5-7.
    
    Parameters:
        T, temperature in K
        P, pressure in Pa
        mWA, molecular weight of component A
        sigmaA, for component A in angstrom
        epsilonA, for component A in K
        mWB, molecular weight of component B
        sigmaB, for component B in angsrom
        epsilonB, for component A in K
        
    Returns:
        Diffusion coefficient of A in B in cm^/s
    """
    k = 1.38064852E-23 # m^2*kg/s^2/K or Pa m^3/K
    NA = 6.022140857E23 # molecules/mol
    #prefix = 3/16*sqrt(4*k*pi*1000*NA)/pi*k*(10**10)**2 100**2 # Pa cm^2/s (g/mol)^(1/2) A^2  
    op = [1.060306,0.15610,0.19300,0.47635,1.03587,1.52996,1.76474,3.89411]
    Ts = k*T/sqrt(epsilonA*epsilonB)
    omegaD = op[0]/(Ts)**op[1]+op[2]*exp(-op[3]*Ts)+op[4]*exp(-op[5]*Ts)+op[6]*exp(-op[7]*Ts)
    mWAB = 2*1/(1/mWA+1/mWB)

    return 266.35232*T**(3/2)/P/mWAB**(1/2)/(1/2*(sigmaA+sigmaB))**2/omegaD

def sigma(cIDs, db):
    """
    Chapman Enskog sigma equations from Perry's Chemical Engineer's Handbook
    Table 5-15, pg 5-49
    
    Parameters:
        cIDs, keys for db
        db, dictionary of database items with keys in cIDs
        
    Returns:
        list of sigma values for each cID in cIDs calculated in various ways...
        item 0 of the list is a two item tuple
        
    """
    sig0 = []
    sig01 = []
    sig1 = []
    sig2 = []
    sig3 = []
    sig4 = []
    for cID in cIDs:
        Tc = db[cID]['critT'] #K
        Pc = db[cID]['critP'] #atm
        Pc = uc(Pc,"bar","atm")
        Vc = db[cID]['critV'] #cm^3/mol
        Tb = db[cID]['nBP'] #K
        Vb = 1/lD(Tb,db[cID]['lDP'])*db[cID]['fW'] #cm^3
        Zc = db[cID]['critC']
        w = db[cID]['acentric']
        Tm = db[cID]['nFP'] #K
        Vm = 1/lD(Tm,db[cID]['lDP'])*db[cID]['fW']
        sig0.append(0.841*Vc**(1/3))
        sig01.append(2.44*(Tc/Pc)**(1/3))
        sig1.append(0.1866*Vc**(1/3)/Zc**(6/5))
        sig2.append(1.18*Vb**(1/3))
        sig3.append(1.222*Vm**(1/3))
        sig4.append((2.3551-0.087*w)*(Tc/Pc)**(1/3))

    return [sig0,sig01,sig1,sig2,sig3,sig4]
        
def epsilonPerk(cIDs, db):
    """
    Table 5-15, pg 5-49 Perry's Chemical Engineer's Handbook
    
    Chapman Enskog epsilon/k equations from Perry's
    """
    epk0 = []
    epk1 = []
    epk2 = []
    epk3 = []
    epk4 = []
    for cID in cIDs:
        Tc = db[cID]['critT']
        Pc = db[cID]['critP']
        Tb = db[cID]['nBP']
        Zc = db[cID]['critC']
        w = db[cID]['acentric']
        Tm = db[cID]['nFP']
        epk0.append(0.77*Tc)
        epk1.append(65.3*Tc*Zc**3.6)
        epk2.append(1.15*Tb)
        epk3.append(1.92*Tm)
        epk4.append((0.7915-0.1693*w)*Tc)

    return [epk0,epk1,epk2,epk3,epk4]

def lJ0(Tc, Vc):
    """
    Lennard Jones sigma and epsilon/k calculation:
        sigma = 0.841*Vc**(1/3)
        epsison/k = 0.75*Tc
    Parameters:
        Tc, critical temperature in K
    Returns:
        [sigma, epsilon/k]
    """
    return [0.841*Vc**(1/3),0.75*Tc]

def lJDB():
    #kludge for __file__ and sys.argv[0] not working when called from exec()
    dbDir = os.path.dirname(os.path.realpath(inspect.getabsfile(lJDB)))
    lJDBPath = os.path.join(dbDir,'..','lennardJones')
    csvFileName = 'lennardJones.csv'
    fileHandle = open(os.path.join(lJDBPath,csvFileName),'r')

    csvReader = csv.DictReader(fileHandle,dialect='excel')
    db = {}
    for row in csvReader:
        db[row['name']] = {}
        for key,val in row.items():
            db[row['name']][key] = val
            
    return db

def DAM(molFrac,DAB):
    """
    multicomponent diffusion coefficient of A in a mix of components m under
    the assumption that the mix of components m all move with the same velocity

    simplification of the Stefan-Maxwell equations
    
    Parameters:
        molFrac, list of mol fractions for each component in DAB excluding
                 component A
        DAB, list of binary diffusion coefficients for A in B where B is a 
             component of mix m
    """
    sumX = 0
    sumXpDAB
    for mF in molFrac:
        sumX += mF
        sumXpDAB += mF/DAB[idx]

    return (1-sumX)/sumXpDAB

if __name__ == "__main__":
    """
    Ref. 1: BSL pg 503, Table 16.2-2
    Ref. 2: "THE PROPERTIES OF GASES AND LIQUIDS", Bruce E. Poling, John M. 
            Prausnitz, John P. O’Connell, 5th ed., McGRAW-HILL, NY, (2001), 
            pgs 11.5-7.

    DAB_data = [[compA, compB, Temp. /K, Pres. /atm, DAB /cm^2/s, ref.]]
    """
    import mpdb
    from mpdb.utils import *
    from mpdb.units import uc, ucTemp
    from mpdb.periodicTable import fW
    from mpdb.pcmpEQ import *
    from tabulate import tabulate

    location = os.path.dirname(os.path.realpath(os.path.abspath(__file__)))
    
    try:
        if not isinstance(db,dict):
            (db, dbMetaData) = dbLoad(interactive=True)
    except:
        (db, dbMetaData) = dbLoad(interactive=True)

    DAB_data = [["carbon dioxide","nitrous oxide",273.2,1.0,0.096,1],
                ["carbon dioxide","carbon monoxide",273.2,1.0,0.139,1],
                ["carbon dioxide","nitrogen",273.2,1.0,0.144,1],
                ["carbon dioxide","nitrogen",288.2,1.0,0.158,1],
                ["carbon dioxide","nitrogen",298.2,1.0,0.165,1],
                ["argon","oxygen",293.2,1.0,0.20,1],
                ["argon","ammonia",333,0.987,0.256,2],
                ["hydrogen","sulfur hexafluoride",298.2,1.0,0.420,1],
                ["hydrogen","methane",298.2,1.0,0.726,1]]
     
    k = 1.38064852E-23 # m^2*kg/s^2/K or Pa m^3/K
    DAB_corr = []
    for idx,dat in enumerate(DAB_data):
        PcA = uc(db[dat[0]]['critP'],dbMetaData[dat[0]]['critP']['units'],"atm")
        TcA = db[dat[0]]['critT']
        mWA = db[dat[0]]['fW']
        PcB = uc(db[dat[1]]['critP'],dbMetaData[dat[1]]['critP']['units'],"atm")
        TcB = db[dat[1]]['critT']
        mWB = db[dat[1]]['fW']
        T = dat[2]
        P = dat[3]
        aDVA = fullerADV(db[dat[0]]['formula'])
        aDVB = fullerADV(db[dat[1]]['formula'])
        DAB0 = birdD(T,P,mWA,PcA,TcA,mWB,PcB,TcB)
        DAB1 = fullerD(T,P*101325,mWA,aDVA,mWB,aDVB)
        DAB_corr.append([DAB0,DAB1])
        
    DABHeader = ['Component A','Component B','Temp. /K','Pres. /atm','Exp. /cm^2/s','Ref.']
    
    print('\n'+tabulate(DAB_data,headers=DABHeader)+'\n')

    DABHeader = ['BSL /cm^2/s','Fuller /cm^2/s']
    print('\n'+tabulate(DAB_corr,headers=DABHeader)+'\n')
    
    ljdb = lJDB()
    cIDs = [key for key in ljdb.keys() if key in db.keys()]
    sig = sigma(cIDs, db)
    eps = epsilonPerk(cIDs,db)
    sig.append([ljdb[cID]['sigma /angstrom'] for cID in cIDs])
    eps.append([ljdb[cID]['epsilonpk /K'] for cID in cIDs])
    
