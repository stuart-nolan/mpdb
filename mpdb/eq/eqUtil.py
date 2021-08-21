"""
mpdb: Material Property Data Base as python module
eqUtil.py: utility methods for thermophysical property estimation

Revision Date: 2021.08.17

SPDX-License-Identifier: BSD-2-Clause
Copyright (c) 2021 Stuart Nolan. All rights reserved.
"""
def func_2arg(arg1,arg2):
    print('arg1 is %s and arg2 is %s' % (arg1,arg2))

def funcDisp(method,*args):
    """
    ref method for function dispatcher from python module
    usage:
      f2arg = funcDisp("2arg")
      f2arg(1,"argument 2") 
      returns:
        arg1 is 1 and arg2 is argument 2
    """
    return getattr(__import__(__name__),'func_%s' % method)

def thermoFuncs(db={},useCoolProp=False):
    """
    add thermo functions from yawsEqn.pyx, coolprop, etc to a db dictionary
    prefer CoolProp functions over yaws pcmpEQ for all CoolProp fluids


    TODO: - CURRENTLY BROKEN while reorganizing "mpdb.eq"
          - fix mCM.py and other demos that depend on "mpdb.eq"
          - sensible strategy to use multiple EQs for same thermo functions
            i.e. CoolProp, pcmpEQ, mcmpEQ for vapor enthalpy...
          - coerce CoolProp into returning value for liquid or vapor states
            given only T like yaws eqns...
    NOTES:
          - python closures are late binding.  The value of variables used in 
            closures are looked up at the time an inner function is called (the 
            lambda's below)
          - python default function arguments are mutable and initialized once
            at the time the function is defined.  Mutalbe means the default 
            argument can be changed at each function call.

    """ 
    cPFluids = coolPropFluids()
    dbtf={}
    for k in db.keys():
        dbtf.update({k: {}})
        if useCoolProp:
            CAS = db[k]['CAS']
            if CAS in cPFluids.keys():
                cPK = cPFluids[CAS]
                dbtf[k].update({'cPKey': cPK})
                dbtf[k].update({'cPEQ': coolPropEqns(cPK)})
                
        if 'vPP' in db[k].keys():
            dbtf[k].update({'vP':lambda T, vPP=db[k]['vPP']: vP(T,vPP)})
        if 'eVP' in db[k].keys():
            dbtf[k].update({'eV':lambda T, eVP=db[k]['eVP']: eV(T,eVP)})
        if 'lDP' in db[k].keys():
            dbtf[k].update({'lD':lambda T, lDP=db[k]['lDP']: lD(T,lDP)})
        if 'lHCP' in db[k].keys():
            dbtf[k].update({'lHC':lambda T, lHCP=db[k]['lHCP']: lHC(T,lHCP)})
        if 'sHCP' in db[k].keys():
            dbtf[k].update({'sHC':lambda T, sHCP=db[k]['sHCP']: sHC(T,sHCP)})
            dbtf[k].update({'sH':lambda T, TRef, sHCP=db[k]['sHCP']: sH(T,sHCP,TRef=TRef)})
        if 'vHCP' in db[k].keys():
            dbtf[k].update({'vHC':lambda T, vHCP=db[k]['vHCP']: vHC(T,vHCP)})
            dbtf[k].update({'vH':lambda T, TRef, vHCP=db[k]['vHCP']: vH(T,vHCP,TRef=TRef)})
        if 'sHCP_f1' in db[k].keys():
            dbtf[k].update({'sHC_f1':lambda T, sHCP=db[k]['sHCP_f1']: hC_f1(T,sHCP)})
            dbtf[k].update({'sH_f1':lambda T, TRef, sHCP=db[k]['sHCP_f1']: h_f1(T,sHCP,TRef=TRef)})
        if 'sHCP_f1p1' in db[k].keys():
            dbtf[k].update({'sHC_f1p1':lambda T, sHCP=db[k]['sHCP_f1p1']: hC_f1(T,sHCP)})
            dbtf[k].update({'sH_f1p1':lambda T, TRef, sHCP=db[k]['sHCP_f1p1']: h_f1(T,sHCP,TRef=TRef)})
        if 'sHCP_f1p2' in db[k].keys():
            dbtf[k].update({'sHC_f1p2':lambda T, sHCP=db[k]['sHCP_f1p2']: hC_f1(T,sHCP)})
            dbtf[k].update({'sH_f1p2':lambda T, TRef, sHCP=db[k]['sHCP_f1p2']: h_f1(T,sHCP,TRef=TRef)})
        if 'vHCP_f1' in db[k].keys():
            dbtf[k].update({'vHC_f1':lambda T, vHCP=db[k]['vHCP_f1']: hC_f1(T,vHCP)})
            dbtf[k].update({'vH_f1':lambda T, TRef, vHCP=db[k]['vHCP_f1']: h_f1(T,vHCP,TRef=TRef)})
        if 'vVP' in db[k].keys():
            dbtf[k].update({'vV':lambda T, vVP=db[k]['vVP']: vV(T,vVP)})
        if 'lVP' in db[k].keys():
            dbtf[k].update({'lV':lambda T, lVP=db[k]['lVP']: lV(T,lVP)})
        if 'vTCP' in db[k].keys():
            dbtf[k].update({'vTC':lambda T, vTCP=db[k]['vTCP']: vTC(T,vTCP)})
        if 'iOTCP' in db[k].keys():
            dbtf[k].update({'iOTC':lambda T, iOTCP=db[k]['iOTCP']: iOTC(T,iOTCP)})
        if 'oLTCP' in db[k].keys():
            dbtf[k].update({'oLTC':lambda T, oLTCP=db[k]['oLTCP']: oLTC(T,oLTCP)})
        if 'sTP' in db[k].keys():
            dbtf[k].update({'sT':lambda T, sTP=db[k]['sTP']: sT(T,sTP)})
        if 'vGEFP' in db[k].keys():
            dbtf[k].update({'vGEF':lambda T,vGEFP=db[k]['vGEFP']: vGEF(T,vGEFP)})
        if 'vEFP' in db[k].keys() and \
           'cPEQ' not in dbtf[k].keys():
            dbtf[k].update({'vEF':lambda T, vEFP=db[k]['vEFP']: vEF(T,vEFP)})
        elif 'eFHUSState' in db[k].keys() and \
             'gas' in db[k]['eFHUSState'] and \
             'cPEQ' in dbtf[k].keys(): \
             dbtf[k].update({'vEF':lambda T, ck=k: db[ck]['eFH298K'] +\
                             dbtf[ck]['cPEQ'].vH(T)/1000})
        elif 'eFHUSState' in db[k].keys() and \
             'gas' in db[k]['eFHUSState'] and \
             'vHCP' in db[k].keys(): \
             dbtf[k].update({'vEF':lambda T, ck=k: db[ck]['eFH298K'] +\
                             dbtf[ck]['vH'](T,298.15)/1000})

    return dbtf
