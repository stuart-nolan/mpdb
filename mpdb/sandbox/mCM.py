"""
mpdb: Material Property Data Base as python module
mCM.py: multi-Component Material

Revision Date: 2021.08.15

SPDX-License-Identifier: BSD-2-Clause
Copyright (c) 2021 Stuart Nolan. All rights reserved.

Notes: this module was previously implemented as a class.  You can go there if 
       you want...
"""
from tabulate import tabulate
import re
import math
import os
from copy import deepcopy
from mpdb.units import uc, ucTemp
from mpdb.utils import dbLoad, dbSubset, dbFind
from mpdb.eq.eqUtil import thermoFuncs

location = os.path.dirname(os.path.realpath(os.path.abspath(__file__)))    

def addVolDens(mCM):
    """
    multi-component density assuming additive volumes

    Parameters:
        mCM, multi-Component Material dictionary - see init

    Returns:
        1/sum([mCM['basisFrac'][cID]/rho for cID, rho in 
               mCM['compDens'].items()])
    """

    if mCM['state'] == 'v':
        mFrac = 'moleFrac'
    else:
        mFrac = 'massFrac'
        
    return 1/sum([mCM[mFrac][cid]/rho for cid, rho in
                  mCM['compDens'].items()])

def conc2compFrac(conc):
    """
    calculate mole or mass fractions from concentrations in mol/vol, mass/vol, 
    or mol/(kg solvent) and sum conc 

    Parameters:
        conc, dictionary like {'<component ID>': concentration} 

    Returns:
        (dict of mole or mass fractions based off conc, sum(conc))
    """
    (concKeys,conc) = compDict2List(conc)
    sumM = sum(conc)
    if sumM > 0:
        compFrac = [c/sumM for c in conc]
    return (compList2Dict(concKeys,compFrac), sumM)

def compDict2List(compDict):
    """
    utility function to create a list of component keys and a list of
    compoenent values from a component dictionary
    """
    return (list(compDict.keys()),list(compDict.values()))

def compList2Dict(compIDs,compVals):
    """
    utility function to create a component dictionary from a list of keys and 
    a list of compoenent values 
    """
    return dict(zip(compIDs,compVals))

def compQ(mCM):
    """
    update mCM['<basis>CompQ'] dictionary where <basis> is mass or mole 

    Parameters:
        mCM, multi-component material dictionary - see init 
    """
    for basis in ['mass', 'mole']:
        basisQ = basis+'Q'
        Q = mCM[basisQ]
        basisCompQ = basis+'CompQ'
        basisFrac = basis+'Frac'
        mCM.update({basisCompQ: dict([(cid,mf*Q) for cid,mf in
                                      mCM[basisFrac].items()])})

def compQ2Frac(mCM):
    """
    update mCM['<basis>Frac'] and mCM['<basis>Q'] from mCM['<basis>CompQ'] and
    then updateBasis(mCM)
    """
    basis = mCM['basis']
    basisCompQ = basis+'CompQ'
    Q = sum(mCM[basisCompQ].values())
    mCM[basis+'Frac'].update(dict([(cid,cQ/Q) for cid, cQ in
                                   mCM[basisCompQ].items()]))
    mCM[basis+'Q'] = Q
    updateBasis(mCM)
    
def dbFTAddVolDens(mCM):
    """
    update mCM['compDens'] and then return addVolDens(mCM), a density 
    assuming additive volumes

    Parameters:
        mCM, multi-Component Material dictionary - see init

    Returns:
        addVolDens(mCM)

    Notes:
        mCM must have the following keys:
        mCM['compDensFunc'] = {'<component>':'<function>', ...}

        mCM['db']['<component>']['<function>'](T) must be a function that
        accepts one parameter (even if that parameter is not used)
    """
    compDens = dict([(cid,mCM['dbtf'][cid][df](mCM['T'])) for cid, df in
                     mCM['compDensFunc'].items()])
    mCM.update({'compDens': compDens})

    return addVolDens(mCM)

def delComp(mCM,delComp=[]):
    """
    delete components from mCM
    """
    delIn = ['moleCompQ','massCompQ', 'moleFrac', 'massFrac', 'conc',
             'compDensFunc', 'compDens']
    for comp in delComp:
        [mCM[item].pop(comp, None) for item in delIn if item in mCM.keys()]
    
def getDens(mCM):
    """
    return density based on type(mCM['dens'])

    mCM['dens'] can be a number or a function

    Parameters:
        mCM, multi-component material dictionary - see init 
    """
    if type(mCM['dens']) in [float, int]:
        dens = mCM['dens']
    else:
        dens = mCM['dens'](mCM)
    return dens

def idGMDens(mCM):
    """
    ideal Gas Molar Density in mol/m^3

    Parameters:
        mCM, multi-component material dictionary - see init 

    Returns
        mCM['P']/R/mCM['T'] in mol/m^3
            
        where

        R = 8.3144599 Pa*m^3/mol/K
    """
    return mCM['P']/8.3144599/mCM['T']
    
def iTM(mCMList=[],ID=''):
    """
    isoThermal Mix (iTM)

    Copy mCMList[0] to oMCM, iTM mCMList[1:] with oMCM and return oMCM.

    If rateUnit is '' in mCMList[0], ignore rateUnit for all items in 
    mCMList.

    If rateUnit is not '' in mCMList[0], use rateUnit for each 
    mCMList[1:] item.  Converting rate units in mCMList[1:] to mCMList[0] rate 
    unit.

    Parameters:
       mCMList, list of 2 or more mCM dicts to combine into one
       ID, string for oMCM['ID']

    Returns:
       oMCM, deepcopy of mCMList[0] + iTM for <basis>Q keys in mCMList[:]

    Assumes the input mCM dicts are "consistent" and "proper":
        - input mCMList <basis>Frac sum to 1
        - all input mCMs have consistent states
        - oMCM will have the same T, P, basis, state, and rateUnit as
          mCMList[0].

    TODO: deepcopy on coolprop functions
    https://stackoverflow.com/questions/7666873/cython-and-deepcopy-woes-with-referenced-methods-functions-any-alternative-id
    """

    if type(mCMList) is list and len(mCMList) > 1:
        for mCM in mCMList:
            if type(mCM) is not dict:
                return "iTM(mCMList): ERROR item (%s) is not dict" % mCM
    else:
        return "iTM(mCMList): ERROR mCMList (%s) is not len(list) > 1" %\
            mCMList
    
    oMCM = deepcopy(mCMList[0])
    oID = oMCM['ID']
    updateBasis(oMCM)
    basis = oMCM['basis']
    basisCompQ = basis+'CompQ'
    basisQUnit = oMCM[basis+'QUnit']
    oCompKeys = list(oMCM[basisCompQ].keys())
    oRateUnit = oMCM['rateUnit']
    if oRateUnit == '': # ignore all rateUnit
        doRate = False
    
    for mCM in mCMList[1:]:
        updateBasis(mCM,toBasis=basis)
        oID = oID + "+" + mCM['ID']
        for cid, cQ in mCM[basisCompQ].items():
            if cid not in oCompKeys:
                oCompKeys.append(cid)
                if doRate:
                    rateUnit = mCM['rateUnit']
                    cQ = uc(cQ,basisQUnit+rateUnit,basisQUnit+oRateUnit)
                    
                oMCM[basisCompQ].update({cid: cQ})
                oMCM[basis+'Frac'].update({cid: mCM[basis+'Frac'][cid]})
            else:
                oMCM[basisCompQ][cid] += cQ

            if cid not in oMCM['db'].keys():
                # oMCM['db'] items linked to mCM['db'] items for now...
                oMCM['db'].update({cid: mCM['db'][cid]})

    if ID == '':
        oMCM['ID']=oID
    else:
        oMCM['ID']=ID

    #create unlinked oMCM['db'] via dbSubset
    oMCM['db'] = dbSubset(oCompKeys,oMCM['db'])
    compQ2Frac(oMCM)

    return oMCM

def init(ID='s0', basis='mass', state='l', T=298.15, dispTUnit='K',
         P=101325, dispPUnit='Pa', moleQ=0, dispMoleQUnit='mole', moleFrac={},
         massQ=0, dispMassQUnit='kg', massFrac={}, rateUnit='',
         aveFw=0, concBasis='molar', conc={},
         dispMolalConcUnit='mol/kg', dispMolarConcUnit='mole/l',
         dispMassConcUnit='kg/l', solvCID='',
         compDensFunc={}, dens=0, densUnit='kg/m^3', dispDensUnit='kg/m^3',
         db={},dbtf={}):
    """
    Initialize a Multi-Component Material (MCM) dictionary
    
    Parameters:
        ID, string ID for mCM  
        state, one of ['s','l','v'] i.e.:
               - (s)olid, (l)iquid, or (v)apor
        basis, one of ['mass', 'mole'].
        T, temperature in K
        dispTUnit, display Unit for T
        P, pressure in Pa
        dispPUnit, display Unit for P
        moleQ, mole quantity in mole
        dispMoleQUnit, display mole unit (do not include /time for rates)
        moleFrac, dict of component mole fractions
        massQ, mass quantity in kg
        dispMassQUnit, display mass unit (do not include /time for rates)
        massFrac, dict of component mass fractions
        rateUnit, rate Unit (/time) i.e. '/s', '/min', '/hr', '/day'
        aveFw, average Formula weight of the mixture in g/mole
        concBasis, one of ['molar', 'mass', 'molal']
        conc, dict of component concentrations
        concUnit, one of ['mole/m^3', 'kg/m^3', 'mole/kg']
        dispMolarConcUnit, mole/vol user input unit for display
        dispMolalConcUnit, mole/kg(solvent) user input unit for display
        dispMassConcUnit, mass/vol user input unit for display
        solvCID, component key that is the solvent for a liquid mixture.
                 'molal' units i.e. 'mol/kg'(solvent).
        compDensFunc, {'<component>':'<function>', ...}.  Function returns a
                      density of the component in densUnit.
        dens, density of mixture/substance.  Default is 0, can be number or
              function (that returns a density) in densUnit.
        densUnit, units for dens.  must be 'mol/m^3' if state = 'v', 
                  'kg/m^3' otherwise
        dispMassDensUnit, display mass density unit (mass/vol)
        dispMoleDensUnit, display mole density unit (mole/vol)
        db, material property database
        dbtf, material property database thermo functions

    Internal computations use the units described below.  db values and 
    functions should be converted to conformable units as necessary.  

    SI base units are used without prefixes:
        kg:        mass;
        m:         length;
        mole:      amount of substance;
        K:         temperature;
        s:         time;
        A:         electric curent;
        cd:        luminous intensity;
    
    Allowed derived SI units or quantities without prefixes:
        /s,
        /min, 
        /hr,
        /day:      rates for mCM['Q'];
        Pa:        Pascal, pressure;
        m^3:       cubic meter, volume;
        J:         Joule, energy;
        C:         columb, electric charge;
        V:         volt, electric potential;
        F:         farad, capacitance;
        ohm:       resistance;

    non SI but used:
        mole/kg(solvent): concentration in mole per kg of solvent;
        mole/m^3   concentration;
        kg/m^3:    concentration;
        g/mole:    formula Weight, aka fW or aveFw;
    """
    mCM = {'ID': ID,
           'T': T,
           'TUnit': 'K',
           'dispTUnit': dispTUnit,
           'P': P,
           'PUnit': 'Pa',
           'dispPUnit': dispPUnit,
           'basis': basis,
           'massQ': massQ,
           'massQUnit': 'kg',
           'dispMassQUnit': dispMassQUnit,
           'massFrac': massFrac,
           'moleQ': moleQ,
           'moleQUnit': 'mole',
           'dispMoleQUnit': dispMoleQUnit,
           'moleFrac': moleFrac,
           'rateUnit': rateUnit,
           'aveFw': aveFw,
           'state': state,
           'conc': conc,
           'dispMolarConcUnit': dispMolarConcUnit,
           'dispMolalConcUnit': dispMolalConcUnit,
           'dispMassConcUnit': dispMassConcUnit,
           'concBasis': concBasis,
           'concUnit': 'mole/m^3', # set by 'concBasis' in sanityCheck
           'solvCID': solvCID,
           'compDensFunc': compDensFunc, 
           'compDens': {},
           'dens': dens,
           'densUnit': densUnit,
           'dispDensUnit': dispDensUnit,
           'db': db,
           'dbtf': dbtf}

    print('\nmCM init: calling sanityCheck(ID: %s)\n---' %
          (mCM['ID'])) 
    sanityCheck(mCM)
    print('---\nmCM init: sanityCheck(ID: %s) done\n'%
          (mCM['ID']))
    
    return mCM
                    
def massFrac2moleFrac(mCM):
    """
    use mCM['massFrac'] to update mCM['moleFrac'] and mCM['aveFw']

    Parameters:
        mCM, multi-component material dictionary - see init 
    """
    aveFw = 0
    db=mCM['db']
    moleFrac = {}
    
    for cid,mf in mCM['massFrac'].items():
        moleFrac.update({cid: mf/db[cid]['fW']})
        aveFw += moleFrac[cid]
        
    if aveFw > 0:
        for cid,mF in moleFrac.items():
            moleFrac[cid] = mF/aveFw
        mCM.update({'aveFw': 1/aveFw})
        mCM.update({'moleFrac': moleFrac})
    else:
        print('massFrac2moleFrac(ID: %s) ERROR: aveFw (%s) =< 0; aborted' %
              (mCM['ID'],aveFw))
        
def sanityCheck(mCM):
    """
    checks for "sane" mCM values
    
    Ref:
    reUnitSep = re.compile(r'(^.*?)(/.*|$)')
    reUnit = reUnitSep.match(mCM['massQUnit']).groups()

    #if mCM['massQUnit'] is 'kg/hr'
    numerator = reUnit[0] # 'kg'
    rateUnit = reUnit[1] # '/hr' 

    if reUnit[0] is not 'kg':
        print("sanityCheck (ID: %s): massQUnit (%s) is not 'kg%s'" %
              (mCM['ID'],mCM['massQUnit'],rateUnit))

    reUnit = reUnitSep.match(mCM['moleQUnit']).groups()        
    if reUnit[0] is not 'mole':
        print("sanityCheck (ID: %s): moleQUnit (%s) is not 'mole'" %
              (mCM['ID'],mCM['moleQUnit']))
    """
    # basis and non zero massQ and/or moleQ
    basis = mCM['basis']
    if basis not in ['mass', 'mole']:
        print("sanityCheck(ID: %s): basis (%s) not 'mass' or 'mole'" %
              (mCM['ID'],basis))
        
    if mCM[basis+'Q'] <= 0:
        print("sanityCheck(ID: %s): %sQ (%s) <= 0" %
              (mCM['ID'],basis,mCM[basis+'Q']))

    # mass and/or mole fraction sum to 1
    tol=1e-9
    sumFrac = sum([mf for mf in mCM['massFrac'].values()])
    if not math.isclose(sumFrac, 1.0, rel_tol=tol):
        print("sanityCheck(ID: %s): sum(massFrac)=%s is not 1 to tol: %s" %
              (mCM['ID'],sumFrac,tol))
        
    sumFrac = sum([mf for mf in mCM['moleFrac'].values()])
    if not math.isclose(sumFrac, 1.0, rel_tol=tol):
        print("sanityCheck(ID: %s): sum(moleFrac)=%s is not 1 to tol: %s" %
              (mCM['ID'],sumFrac,tol))

    # density unit conforms to mCM convention
    if type(mCM['dens']) not in [float, int]: # if not a constant density
        if mCM['state'] == 'v':
            # force density unit to mole/m^3 if state = 'v'
            mCM['densUnit'] = 'mole/m^3'
        else:
            # force density unit to kg/m^3 for everything else
            mCM['densUnit'] = 'kg/m^3'
    else: # constant density, warn if dens<=0 and/or densUnit non-conforming
        if mCM['dens'] <= 0:
                print("sanityCheck(ID: %s): mCM['dens']=%s <= 0" %
                      (mCM['ID'],mCM['dens']))            

        if mCM['state'] == 'v':
            # densUnit should be mole/m^3 if state = 'v'
            if mCM['densUnit'] != 'mole/m^3':
                print("sanityCheck(ID: %s): state='v', mCM['densUnit'] (%s) is not 'mole/m^3'" % (mCM['ID'],mCM['densUnit']))
        else:
            # density unit should be kg/m^3 for everything else
            if mCM['densUnit'] != 'kg/m^3':
                print("sanityCheck(ID: %s): state='%s', mCM['densUnit'] (%s) is not 'kg/m^3'" % (mCM['ID'],mCM['state'],mCM['densUnit']))

    if mCM['state'] == 'v': # check dispDensUnit
        reUnitSep = re.compile(r'(^.*?)(/.*|$)')
        reUnit = reUnitSep.match(mCM['dispDensUnit']).groups()
        if reUnit[0] not in 'mole':
            print("sanityCheck(ID: %s): state='v', mCM['dispDensUnit'] (%s) is not 'mole%s'" % (mCM['ID'],mCM['dispDensUnit'],reUnit[1]))

    # concBasis is 'molal', 'molar', or 'mass' and set concUnit accordingly
    if mCM['concBasis'] == 'molal':
        mCM.update({'concUnit': 'mole/kg'})
    elif mCM['concBasis'] == 'mass':
        mCM.update({'concUnit': 'kg/m^3'})
    elif mCM['concBasis'] == 'molar':
        mCM.update({'concUnit': 'mole/m^3'})
    else:
        print("sanityCheck(ID: %s): mCM['concBasis'] (%s) is not 'molal', 'mass' or 'molal'" % (mCM['ID'],mCM['concBasis']))


def moleFrac2massFrac(mCM):
    """
    use mCM['moleFrac'] to update mCM['massFrac'] and mCM['aveFw']

    Parameters:
        mCM, multi-component material dictionary - see init 
    """
    
    aveFw = 0
    db=mCM['db']
    massFrac = {}
    
    for cid, mf in mCM['moleFrac'].items():
        massFrac.update({cid: mf*db[cid]['fW']})
        aveFw += massFrac[cid]
        
    if aveFw > 0:
        for cid, mF in massFrac.items():
            massFrac[cid] = mF/aveFw
        mCM.update({'massFrac': massFrac})
        mCM.update({'aveFw': aveFw})
    else:
        print('moleFrac2massFrac(ID: %s) ERROR: aveFw (%s) =< 0; aborted' %
              (mCM['ID'],aveFw))
            
def pTable(mCM, tables=['prop','comp','conc'], dispUnit=True, updateBasis=False):
    """
    print mCM info in tables
    
    Parameters:
        mCM, multi-component material dictionary - see init 
        tables, list or string containing ['comp', 'prop', 'conc']
        dispUnits, True to use display units or False to use internal units
    """    
    basis = mCM['basis']
    Basis = mCM['basis'].capitalize()
    if updateBasis:
        updateBasis(mCM)
    elif basis == 'mole':
        moleFrac2massFrac(mCM)
    elif basis == 'mass':
        massFrac2moleFrac(mCM)
    
    aveFw = mCM['aveFw']

    if dispUnit:
        updateDispUnit(mCM)
        T = mCM['dispT']
        TUnit = mCM['dispTUnit']
        P = mCM['dispP']
        PUnit = mCM['dispPUnit']
        Q = mCM['disp'+Basis+'Q']
        QUnit = mCM['disp'+Basis+'QUnit']+mCM['rateUnit']
        conc = mCM['dispConc']
        concUnit = mCM['disp'+mCM['concBasis'].capitalize()+'ConcUnit']
        dens = mCM['dispDens']
        densUnit = mCM['dispDensUnit']
    else:
        T = mCM['T']
        TUnit = mCM['TUnit']
        P = mCM['P']
        PUnit = mCM['PUnit']
        Q = mCM[basis+'Q']
        QUnit = mCM[basis+'QUnit']
        conc = mCM['conc']
        concUnit = mCM['concUnit']
        dens = getDens(mCM)
        densUnit = mCM['densUnit']

    showVol = False
    if dens > 0:
        reUnitSep = re.compile(r'(^.*/)(.*|$)')
        reUnit = reUnitSep.match(densUnit).groups()

        vol = Q/dens
        showVol = True
        if dispUnit:
            volUnit = reUnit[1]+mCM['rateUnit']
        else:
            volUnit = reUnit[1]
    
    if 'prop' in tables:
        propHeader = ['%s Property' % mCM['ID'],'Value','Units']
        propTable =  [['Basis',basis,'-'],
                      ['State',mCM['state'],'-'],
                      ['T',"{0:.4E}".format(T),TUnit],
                      ['P',"{0:.4E}".format(P),PUnit],
                      ['Quantity',"{0:.4E}".format(Q), QUnit],
                      ['Density',"{0:.4E}".format(dens), densUnit],
                      ['Ave Fw',"{0:.4E}".format(aveFw),'g/mole']]
        if showVol:
            propTable.insert(5,['Volume',"{0:.4E}".format(vol), volUnit])
            
        print('\n'+tabulate(propTable,headers=propHeader)+'\n')
        
    if 'comp' in tables:
        compHeader = ['%s Component ID' % mCM['ID'],
                      'Mass Fraction',
                      'Mole Fraction',
                      'Quantity /%s' % QUnit]        
        compTable = []
        for cid, mf in mCM[basis+'Frac'].items():
            compTable.append([cid,
                              "{0:.4E}".format(mCM['massFrac'][cid]),
                              "{0:.4E}".format(mCM['moleFrac'][cid]),
                              "{0:.4E}".format(mf*Q)])
        print('\n'+tabulate(compTable,headers=compHeader)+'\n')

    if 'conc' in tables:
        concHeader = ['%s Component ID' % mCM['ID'],
                      'Mass Fraction',
                      'Mole Fraction',
                      'Conc /%s' % concUnit]
        concTable = []
        for cid, cval in conc.items():
            concTable.append([cid,
                              "{:.4E}".format(mCM['massFrac'][cid]),
                              "{:.4E}".format(mCM['moleFrac'][cid]),
                              "{:.4E}".format(cval)])
        print('\n'+tabulate(concTable,headers=concHeader)+'\n')

def updateBasis(mCM,toBasis=''):
    """
    update mCM['basis'], mCM['massFrac'], mCM['moleFrac'], mCM['massQ'],
    mCM['moleQ'], and then call compQ(mCM) and updateConc(mCM)

    Parameters:
        mCM, multi-component material dictionary - see init
        toBasis, set mCM['basis']=toBais if toBasis in ['mole', 'mass'] 
    """
    basis = mCM['basis']
        
    if basis == 'mass':
        massFrac2moleFrac(mCM)
        mCM['moleQ'] = mCM['massQ']*1000/mCM['aveFw']
        
    if basis == 'mole':
        moleFrac2massFrac(mCM)
        mCM['massQ'] = mCM['moleQ']*mCM['aveFw']/1000

    compQ(mCM)
    # TODO: should updateConc be called here?
    updateConc(mCM)
    
    if toBasis in ['mass', 'mole']:
        mCM['basis'] = toBasis        

def updateByCompQ(mCM):
    """
    """
    basis = mCM['basis']
    basisCompQ = basis+'CompQ'
    basisFrac = basis+'Frac'
    basisQ = basis+'Q' 
    mCM[basisQ] = sum([comp for comp in mCM[basisCompQ].values()]) 
    mCM[basisFrac] = dict([(cid, compQ/mCM[basisQ]) for cid,compQ in
                           mCM[basisCompQ].items()])
    updateBasis(mCM)
    
def updateConc(mCM):
    """
    update mCM['conc'] 

    Calculate component concentrations from mass or mol based quantities and
    the density of the component mixture.  Density from mCM['dens'].

    set mCM['concBasis'] before calling to set/change conc basis.

    Parameters:
        mCM, multi-component material dictionary - see init    
    """    
    rho = getDens(mCM)
    if rho <= 0:
        print("updateConc (ID: %s): mCM['dens']=%s <= 0; aborted" %
              (mCM['ID'],mCM['dens']))
        return

    if mCM['state'] == 'v': 
        # mCM convension is density in mol/m^3 if state is 'v'
        # but updateConc requires rho in kg/m^3
        rho = rho*mCM['aveFw']/1000 # kg/m^3
        
    if mCM['concBasis'] == 'molal':
        # need mole component and kg solvent
        mCM['concUnit'] = 'mole/kg'
        cQ = mCM['moleCompQ']
        if mCM['basis'] == 'mole':
            # kg solvent
            invDenom = 1000*mCM['db'][mCM['solvCID']]['fW']/ \
                       cQ[mCM['solvCID']]
        elif mCM['basis'] == 'mass':
            # kg solvent
            invDenom = 1/(mCM['massQ']*mCM['massFrac'][mCM['solvCID']])
    elif mCM['concBasis'] == 'mass':
        # need kg component and density
        mCM['concUnit'] = 'kg/m^3'
        cQ = mCM['massCompQ']
        if mCM['basis'] == 'mole':
            invDenom = rho/mCM['moleQ']*mCM['aveFw']/1000
        elif mCM['basis'] == 'mass':
            invDenom = rho/mCM['massQ']
    elif mCM['concBasis'] == 'molar':
        # need mole component, aveFw, and density*1000
        mCM['concUnit'] = 'mole/m^3'
        cQ = mCM['moleCompQ']
        if mCM['basis'] == 'mole':
            invDenom = rho*1000/(mCM['moleQ']*mCM['aveFw'])
        elif mCM['basis'] == 'mass':
            invDenom = rho/mCM['massQ']
        
    mCM.update({'conc': dict([(cid,cQval*invDenom)
                                  for cid, cQval in cQ.items()])})

def updateDispUnit(mCM):
    """
    update mCM['dispT'], mCM['dispP'], mCM['dispBasisQ'], mCM['dispDens'], and
    mCM['dispConc']
    """
    mCM.update({'dispT': ucTemp(mCM['T'],
                                mCM['TUnit'],
                                mCM['dispTUnit'])})
    mCM.update({'dispP': uc(mCM['P'],
                            mCM['PUnit'],
                            mCM['dispPUnit'])})

    for basis in ['mole','mass']:
        Basis = basis.capitalize()
        Q = mCM[basis+'Q']
        fromQUnit = mCM[basis+'QUnit']
        toQUnit = mCM['disp'+Basis+'QUnit']
        mCM.update({'disp'+Basis+'Q': uc(Q,fromQUnit,toQUnit)})

    mCM.update({'dispDens': uc(getDens(mCM),
                               mCM['densUnit'],
                               mCM['dispDensUnit'])})

    toConcUnit = mCM['disp' + mCM['concBasis'].capitalize() + 'ConcUnit']
    mCM.update({'dispConc':
                dict([(cid,
                       uc(conc,
                          mCM['concUnit'],
                          toConcUnit))
                      for cid, conc in mCM['conc'].items()])})

if __name__ == "__main__":
    """
    """
    try:
        if not isinstance(db,dict):
            (db, dbMetaData) = dbLoad(interactive=True)
    except:
        (db, dbMetaData) = dbLoad(interactive=True)

    s2cID = ['n-hexane','n-heptane']
    s2massF = [0.4, 0.6]
    s2db=dbSubset(s2cID,db)
    
    # define how density (or other thermo funcs are calculated)
    s2dbtf=thermoFuncs(db=s2db)
    s2compDensFunc = dict([(cID, 'lD') for cID in s2cID if 'lD' in
                           s2dbtf[cID].keys()])
    
    s2=init(ID='s2',
            state='l',
            T=288.15, # K
            P=101325, # Pa
            massFrac=compList2Dict(s2cID,s2massF),
            basis='mass',
            massQ=100, # kg
            dispTUnit='F',
            dispPUnit='psi',
            dispMassQUnit='lb',
            rateUnit='/min',
            dispDensUnit='lb/gal',
            compDensFunc=s2compDensFunc,
            dens=dbFTAddVolDens, # comment out to test dens=0
            #densUnit='g/ml', # test non-conforming dens unit
            concBasis='mass',
            dispMassConcUnit='lb/gal',
            db=s2db,
            dbtf=s2dbtf)

    updateBasis(s2)
    pTable(s2)
