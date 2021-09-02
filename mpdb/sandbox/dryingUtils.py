"""
SPDX-License-Identifier: BSD-2-Clause
Copyright (c) 2021 Stuart Nolan. All rights reserved.
"""
from math import log
import os
from mpdb.eq import yaws_python as yeq
from mpdb.eq import misc_python as meq
from mpdb.periodicTable import fW, elemCount
import mpdb.sandbox.diffusionCoef as dC
import mpdb.sandbox.mCM as mCM
from mpdb.units import uc, ucTemp
from mpdb.utils import dbLoad, dbSubset, dbFind

def Le(tC,hC,rho,DAB):
    """
    Lewis number: 
        thermal conductivity/(heat capacity)/(humid vapor density)/(DAB)
    
    Paramters:
        tC, thermal conductivity in W/m/K
        hC, heat capacity in J/mol/K
        rho, molar density of humid vapor in mol/m^3
        DAB, diffusion of component A in B in m^2/s

    Returns:
        Lewis number (dimensionless)
        
    """
    return tC/hC/rho/DAB
    
def combustRxnStoic(fuelKeys,db={}):
    """
    combustion Reaction Stoichiometry

    for a fuel reacting with oxygen according to the reaction:

        CxHyOzSw + [(4x + y - 2z + 4w)/4]O2 -> xCO2 + (y/2)H2O + wSO2

    Parameters:
        fuelKeys, list of component keys in db that react with oxygen 
        according to the reaction above.  These components must have a 
        formula consistent with CxHyOzSw and the defined reaction.

        db, material property database dict containing all keys in fuelKeys

    Returns:
        rxns, dict of reactions with nested 'reactants' and 'products' 
              dictionaries for each key in fuelKeys.
    
    Example:
        (db, dbMetaData) = dbLoad(interactive=True)
        natGas = ['methane','ethane','propane','butane','hydrogen',
                  'hydrogen sulfide']
        nGRS = combustRxnStoic(natGas,db=db)
        nGRS
        {'methane': {'reactants': {'methane': 1,
                                   'oxygen': 2},
                     'products': {'carbon dioxide': 1,
                                  'water': 2}},
         'ethane': {'reactants': {'ethane': 1,
                                  'oxygen': 14/4},
                    'products': {'carbon dioxide': 2,
                                 'water': 3}},
         'propane': {'reactants': {'propane': 1,
                                   'oxygen': 5},
                     'products': {'carbon dioxide': 3,
                                  'water': 4}},
         'hydrogen': {'reactants': {'hydrogen': 1,
                                    'oxygen': 1/2},
                      'products': {'water': 1/2}},
         'hydrogen sulfide': {'reactants': {'hydrogen sulfide': 1,
                                            'oxygen': 2},
                              'products': {'sulfur dioxide': 1,
                                           'water': 1}}}
    """
    ignoreCIDs = ['oxygen','carbon dioxide','sulfur dioxide','water']
    rxns = {}
    for cid in fuelKeys:
        cidFormula = db[cid]['formula']
        if cid in ignoreCIDs:
            print("combustRxns: key '%s' in fuelKeys with formula '%s' ignored" % (cid,cidFormula))
            continue
        eC = elemCount(cidFormula)
        # simple filter for components that have elements not accepted
        if len([key for key in eC.keys() if key not in
               ['C', 'H', 'O', 'S']]) != 0:
            print("combustRxns: key '%s' in fuelKeys with formula '%s' ignored" % (cid,cidFormula))
            continue
        prodDict = {}
        x = y = z = w = 0
        if 'C' in eC:
            x = eC['C']
            prodDict.update({'carbon dioxide': x})
        if 'H' in eC:
            y = eC['H']
            prodDict.update({'water': y/2})
        if 'O' in eC:
            z = eC['O']
        if 'S' in eC:
            w = eC['S']
            prodDict.update({'sulfur dioxide': w})

        rxns.update({cid: {'reactants': {cid: 1,
                                         'oxygen': (4*x+y-2*z+4*w)/4},
                           'products': prodDict}})        
    return rxns
    
def combustNaturalGas(nG,
                      air=mCM.init(ID='air:O2+N2',
                                   state='v',
                                   basis='mole',
                                   moleQ=1,
                                   moleFrac={'oxygen': 0.21,
                                             'nitrogen': 0.79},
                                   dens=mCM.idGMDens,
                                   dispDensUnit='mole/l'),
                      excessAir=0,
                      nGRS={}):
    """
    Combust natural gas with air according to the reaction:

        CxHyOzSw + [(4x + y - 2z + 4w)/4]O2 -> xCO2 + (y/2)H2O + wSO2

    Parameters:
        nG, multi-component material dict (see initMCM) containing at least 
            one component with a formula consistent with CxHyOzSw.  Note
            fuelKeys below.  Combustion reaction is assumed to go to
            completion.

        air, multi-component material dictionary containing an O2 component,
             see initMCM.  Defaults is O2: .21 and N2: 0.79.  air['moleQ'],
             air['massQ'], air['moleCompQ'], air['massCompQ'], and 
             air['rateUnit'] will be changed based on nG and excessAir.
             
        excessAir, percent (0-100) excess air requested over theoretical air 
                   needed to completly combust nG.

        nGRS, dict - natural gas combustion Stoichiometry. See combustRxnStoic.

    Return:
        cProd, combustion products in a (new) mCM dict

    """
    mCM.updateBasis(nG,toBasis='mole')
    mCM.updateBasis(air,toBasis='mole')
    
    rxnCompQ = {'oxygen': 0,
                'carbon dioxide': 0,
                'water': 0,
                'sulfur dioxide': 0}
    nGDelKey = []
    for cid, Q in nG['moleCompQ'].items():
        if cid in nGRS.keys():
            # amount reactants consumed and products produced per mole nG feed
            rxnCompQ['oxygen'] += Q*nGRS[cid]['reactants']['oxygen']
            for key in ['carbon dioxide','water','sulfur dioxide']:
                if key in nGRS[cid]['products'].keys():
                    rxnCompQ[key] += Q*nGRS[cid]['products'][key]
            # list of reactants to delete in cProd as the reaction is assumed
            # to go to completion
            nGDelKey.append(cid)

    # modify air, iTM air with nG, and then update cProd
    air['moleQ'] = rxnCompQ['oxygen']*(1+excessAir/100)/\
                   air['moleFrac']['oxygen']
    mCM.updateBasis(air)
    cProd=mCM.iTM([nG,air],ID='cNG')
    # remove fuelKeys and oxygen from cProd if excessAir = 0
    mCM.delComp(cProd,nGDelKey)
    if excessAir == 0:
        mCM.delComp(cProd,['oxygen'])
    else:
        cProd['moleCompQ']['oxygen'] -= rxnCompQ['oxygen']
        
    for cid in ['carbon dioxide','water','sulfur dioxide']:
        if cid in cProd['moleCompQ'].keys():
            cProd['moleCompQ'][cid] += rxnCompQ[cid]

    mCM.updateByCompQ(cProd)
           
    return cProd

if __name__ == "__main__":
    try:
        if not isinstance(db,dict):
            (db, dbmd) = dbLoad(interactive=True)
    except:
        (db, dbmd) = dbLoad(interactive=True)
    # dry air mole fractions
    dAMoleFrac = {'nitrogen': 0.78084,
                  'oxygen': 0.20946,
                  'argon': 0.0093,
                  'carbon dioxide': 0.0004}
    
    evap = 'water'
    yEvap = 0.01 # mole fraction of evap comp, water RH ~ 50% @ 298.15 K
    yEvapM1 = 1-yEvap

    # humid air mole fractions
    hAMoleFrac = dict([(cid, yEvapM1*mf) for cid, mf in
                       dAMoleFrac.items()])
    hAMoleFrac.update({evap: yEvap})
    (hACIDList,hAMoleFracList) = mCM.compDict2List(hAMoleFrac)
    hAdb = dbSubset(hACIDList,db)
    
    hA = mCM.init(ID='hA',
                  basis='mole',
                  state='v',
                  T=298.15, # K
                  dispTUnit='F',
                  P=101325, # Pa
                  dispPUnit='psia',
                  moleQ=1,
                  moleFrac=hAMoleFrac,
                  dens=mCM.idGMDens,
                  dispDensUnit='mole/l',
                  db=hAdb)
    mCM.updateBasis(hA)
    hA.update({'dbtf': meq.thermoFuncs(db=hAdb)})
    mCM.pTable(hA)

    vVL = [yeq.vaporViscosity(hA['T'],db[cID]['vVP']) for cID in hACIDList]
    vTCL = [yeq.vaporThermalConductivity(hA['T'],db[cID]['vTCP'])
            for cID in hACIDList]
    vHCL = [yeq.vaporHeatCapacity(hA['T'],db[cID]['vHCP'])
            for cID in hACIDList]
    fWL = [db[cID]['fW'] for cID in hACIDList]
    hATC = meq.mCVTC(hAMoleFracList,vTCL,vVL,fWL)
    hAHC = meq.mFWQ(hAMoleFracList,vHCL)
    evapADV = dC.fullerADV(db[evap]['formula'],ring=False)
    airADV = dC.fullerADV('air')
    fWA = db[evap]['fW']
    fWB = 28.97
    DAB = dC.fullerD(hA['T'],hA['P'],fWA,evapADV,fWB,airADV)
    DAB = uc(DAB,"cm^2/s","m^2/s")
    rhohA = hA['P']/8.314/hA['T'] # mol/m^3
    hALe = Le(hATC,hAHC,rhohA,DAB)
    Y = hA['massFrac'][evap]/(1-hA['massFrac'][evap])
    Ysat = uc(yeq.vaporPressure(hA['T'],db[evap]['vPP']),'mmHg','Pa')/hA['P']
    Ysat = Ysat*fWA/((1-Ysat)*fWB)
    phi = fWA/fWB/(Ysat-Y)*log(1+(Ysat-Y)/(fWA/fWB+Y))

    from CoolProp.CoolProp import PropsSI as pSI
    cPMF = [yEvap,yEvapM1]
    cPT = hA['T'] #K
    cPP = hA['P'] #Pa
    fWA = pSI('M',evap)*1000
    fWB = pSI('M','Air')*1000
    cPTC = [pSI('L','T',cPT,'Q',1,evap),pSI('L','T',cPT,'P',cPP,'Air')]
    cPHC = [pSI('Cpmolar','T',cPT,'Q',1,evap),
            pSI('Cpmolar','T',cPT,'P',cPP,'Air')]
    cPV = [pSI('V','T',cPT,'Q',1,evap),pSI('V','T',cPT,'P',cPP,'Air')]
    cPHCm = yEvap*cPHC[0]+yEvapM1*cPHC[1]
    cPTCm = meq.mCVTC(cPMF,cPTC,cPV,[fWA,fWB])
    cPDAB = dC.fullerD(cPT,cPP,fWA,evapADV,fWB,airADV)
    cPDAB = uc(cPDAB,"cm^2/s","m^2/s")
    cPRho = cPP/8.314/cPT # mol/m^3
    cPLe = Le(cPTCm,cPHCm,cPRho,cPDAB)
    cPY = yEvap*fWA/(yEvapM1*fWB)
    cPYsat = pSI('P','T',cPT,'Q',0,evap)/cPP
    cPYsat = cPYsat*fWA/((1-cPYsat)*fWB)
    cPphi = fWA/fWB/(cPYsat-cPY)*log(1+(cPYsat-cPY)/(fWA/fWB+cPY))
    
    result=(hALe**(-2./3)*phi,cPLe**(-2./3)*cPphi)

    print("Lewis number for %s in air at %s %s and %s %s:" % (evap,
                                                              hA['T'],
                                                              hA['TUnit'],
                                                              hA['P'],
                                                              hA['PUnit']))
    print("yaws: %s; coolprop: %s;\n\n" % result)

        
    print("Combust natural gas in humid air\n")
    """
    natGas, assumed near the midpoint of perry's table 27.8 "Analysis of 
            Natural Gas" as a dictionary:
    """
    
    nGMoleFrac={'methane': 0.9095,
                'ethane': 0.053,
                'propane': 0.017,
                'nitrogen': 0.015,
                'carbon dioxide': 0.0055}
    nGdb = dbSubset(mCM.compDict2List(nGMoleFrac)[0],db)
                        
    nGComps = list(nGMoleFrac.keys())
    nGRS = combustRxnStoic(nGComps,nGdb)
    
    natGas = mCM.init(ID='natGas',
                      state='v',
                      basis='mole',
                      T=298.15, # K
                      dispTUnit='F',
                      P=101325, # Pa
                      dispPUnit='psia',
                      moleQ=1,
                      moleFrac=nGMoleFrac,
                      dens=mCM.idGMDens,
                      dispDensUnit='mole/l',
                      db=nGdb)

    mCM.updateBasis(natGas)
    natGas.update({'dbtf': meq.thermoFuncs(db=nGdb)})
    mCM.pTable(natGas)
    
    """
    Nist Shomate equation implementation (in mpdb/eq/meq.pyx) 
    validation with co2
    REF: https://webbook.nist.gov/cgi/cbook.cgi?ID=C124389&Units=SI&Mask=1&Type=JANAFG&Table=on#JANAFG
    """
    CO2_shomate_TBounds=[[298.0, 1200], [1200, 6000]] # K
    CO2_shomate_params=[[24.99735, 55.18696, -33.69137, 7.948387, 
                         -0.136638, -403.6075, 228.2431, -393.5224],
                        [58.16639, 2.720074, -0.492289, 0.038844, -6.447293,
                         -425.9186, 263.6125, -393.5224]]
    db['carbon dioxide'].update({'sEQP': CO2_shomate_params})
    db['carbon dioxide'].update({'sEQTB': CO2_shomate_TBounds})
    
    H2O_shomate_TBounds=[[500, 1700], [1700, 6000]] # K
    H2O_shomate_params=[[30.09200, 6.832514, 6.793435, -2.534480,
                         0.082139, -250.8810, 223.3967, -241.8264],
                        [41.96426, 8.622053, -1.499780, 0.098119,
                         -11.15764, -272.1797, 219.7809, -241.8264]]
    db['water'].update({'sEQP': H2O_shomate_params})
    db['water'].update({'sEQTB': H2O_shomate_TBounds})

    N2_shomate_TBounds = [[100., 500.], [500., 2000.], [2000., 6000.]]
    N2_shomate_params = [[28.98641, 1.853978, -9.647459, 16.63537,
                          0.000117, -8.671914, 226.4168, 0.0],
                         [19.50583, 19.88705, -8.598535, 1.369784,
                          0.527601, -4.935202, 212.39, 0.0],
                         [35.51872, 1.128728, -0.196103, 0.014662,
                          -4.55376, -18.97091, 224.981, 0.0]]
    db['nitrogen'].update({'sEQP': N2_shomate_params})
    db['nitrogen'].update({'sEQTB': N2_shomate_TBounds})

    Ar_shomate_TBounds = [[298., 6000.]] # K
    Ar_shomate_params = [[20.78600, 2.825911E-07, -1.464191E-07,
                          1.092131E-8, -3.661371E-8, -6.197350,
                          179.9990, 0.000000]]
    db['argon'].update({'sEQP': Ar_shomate_params})
    db['argon'].update({'sEQTB': Ar_shomate_TBounds})

    CO2_validation_data=[["Temperature /K","Cp /J/mol/K","S /J/mol/K",
                          "H/kJ/mol"],
                         [298.,37.12,213.8,-0.01], [300.,37.22,214.0,0.07],
                         [400.,41.34,225.3,4.00], [500.,44.61,234.9,8.31],
                         [600.,47.32,243.3,12.91], [700.,49.57,250.8,17.75],
                         [800.,51.44,257.5,22.81], [900.,53.00,263.6,28.03],
                         [1000.,54.30,269.3,33.40],
                         [1100.,55.40,274.5,239.2,38.89],
                         [1200.,56.35,279.4,242.3,44.47]]
    
    CO2_shomate_enthalpy=[meq.enthalpy_Shomate_bounded(val[0],
                                                       CO2_shomate_TBounds,
                                                       CO2_shomate_params)
                          for val in CO2_validation_data[1:]]

    cNG = combustNaturalGas(natGas,air=hA,excessAir=0,nGRS=nGRS)
    cNG.update({'dbtf': meq.thermoFuncs(db=cNG['db'])})
    #cNG.update({'dbtf': meq.thermoFuncs(db=cNG['db'],useCoolProp=True)})
    #mCM.pTable(cNG)
    
    TRef = natGas['T']
    for cid, cQ in cNG['moleCompQ'].items():
        cNG['dbtf'][cid].update({'vH':lambda T, TRef, sEQTB=db[cid]['sEQTB'], sEQP=db[cid]['sEQP']: meq.enthalpy_Shomate_bounded(T,sEQTB,sEQP,TRef=TRef)})
    # Energy balance via vapor enthalpy of formation (vEF)
    # TRef chosen from nG input stream.
    #natGas.update({'dbtf': meq.thermoFuncs(db=natGas['db'],useCoolProp=True)})
    #hA.update({'dbtf': meq.thermoFuncs(db=hA['db'],useCoolProp=True)})    
    hAEF = sum([cQ*hA['dbtf'][cid]['vEF'](TRef) for cid,cQ in
                hA['moleCompQ'].items()])
    nGEF = sum([cQ*natGas['dbtf'][cid]['vEF'](TRef) for cid,cQ in
                natGas['moleCompQ'].items()])
    cNGEF = lambda T, P0=cNG,P1=hAEF,P2=nGEF:\
        sum([cQ*(db[cid]['sEQP'][0][7] +
                 P0['dbtf'][cid]['vH'](T,TRef))
             for cid, cQ in
             P0['moleCompQ'].items()])-P1-P2
    
    from scipy.optimize import fsolve
    (T,infodict,ler,mesg)=fsolve(cNGEF,1900,full_output=1,factor=1)
    cNG['T'] = float(T)
    mCM.pTable(cNG)
    print("Temperature of combustion gas is a calculated adiabatic flame\ntemperature (dissociation not considered).  \nMeasured value: 1960 deg. C, 3562 deg. F\nTheoretical (dissociation not considered) methane in dry air \nat 1 atm: 2054 deg. C, 3729 deg. F\nRef: http://elearning.cerfacs.fr/combustion/n7masterCourses/adiabaticflametemperature/index.php")

    cNGE=sum([mf*cNG['dbtf'][cid]['vH'](cNG['T'],TRef) for
              cid, mf in cNG['moleFrac'].items()])
