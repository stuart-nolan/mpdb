"""
mpdb: Material Property Data Base as python module
rSanderImport.py: import Henry's law data compiled by R. Sander

Henry's law data from R. Sander, "Compilation of Henry’s law constants
(version 4.0) for water as solvent", Atmos. Chem. Phys., 15,
4399–4981, 2015

Revision Date: 2021.08.15

SPDX-License-Identifier: BSD-2-Clause
Copyright (c) 2021 Stuart Nolan. All rights reserved.
"""

import os
import re
import csv
from time import strftime
import mpdb

def rSander2mpdb(rSdb,db,dbmd):
    """
    Parameters:
      rSdb: dictionary returned by importHbpSIf90(cfg)
      db: mpdb compatible dictionary
      dbmd: mpdb compatible dictionary of meta data
    """
    from mpdb.pcmpEQ import hbpSI
    from mpdb.periodicTable import fW
    refndSpace = re.compile(r'\s+')

    importDateTime = strftime("%Y%m%d-%H%M")
    for key in rSdb.keys():
        if key not in db.keys():
            db.update({key:{'name':key,
                              'CAS':rSdb[key]['CAS'],
                              'formula':rSdb[key]['formula'],
                              'fW':fW(rSdb[key]['formula'])}})
        db[key].update({'hbpSIP':[],
                          'hbpSIPL':rSdb[key]['henryData'],
                          'hbpSI_index':None})
        dbmd.update({key:{'hbpSI':
                                  {'units':'mol/kg/Pa',
                                   'description':'Henry\'s law data & equation from R. Sander',
                                   'eqnDoc':hbpSI.__doc__,
                                   'dataSource':cfg['hbpSIf90FileName'],
                                   'importDateTime':importDateTime}}})

def importHbpSIf90(cfg):
    """ 
    Load Henry's law equation data from the rsander/data/HbpSI.f90 file.
        
    HbpSI.f90 is in molality per pressure units, 
    i.e. mol solute/(kg water * Pa)
    """
    reWS = re.compile(r'\s+')
    reSpecies = re.compile(r'! species:\s+(.*?)$')
    reCAS = re.compile(r'! casrn:\s+(.*?)?$')
    reFormula = re.compile(r'! formula:\s+(.*?)$')
    reFPItems = re.compile(r'(.*?)(\(.*?\)_[0-9]*)')
    reFPEl = re.compile(r'\((.*?)\)_([0-9]*)')
    reFEls = re.compile(r'([A-Z][a-z]?)_*\{*([0-9]*)\}*')
    reHD = re.compile(r'.*? = *(.*?)\s+\D*([0-9]*).*?!.*?type: (.*?), ref: (.*)$')
    with open(cfg['hbpSIf90FQPN'],'r') as hDF:
        db = {}
        for line in hDF:
            if reHD.match(line) is not None:
                db[specKey]['henryData'].append(reHD.match(line).groups())
            elif reSpecies.match(line) is not None:
                specKey = reSpecies.match(line).groups()[0]
                db.update({specKey:{}})
                db[specKey].update({'henryData':[]})
            elif reCAS.match(line) is not None:
                db[specKey].update({'CAS':reCAS.match(line).groups()[0]})
            elif reFormula.match(line) is not None:
                formula = reFormula.match(line).groups()[0]
                fList = []
                if reFPItems.findall(formula):
                    for pItem in reFPItems.findall(formula):
                        if pItem[0] != '':
                            fList.extend([[item[0],int(item[1])] if item[1] != '' else [item[0],1] for item in reFEls.findall(pItem[0])])
                        fPEl = reFPEl.match(pItem[1]).groups()
                        fEls = [list(item) for item in reFEls.findall(fPEl[0])]
                        for idx, fEl in enumerate(fEls):
                            if fEl[1] == '': 
                                fEl[1] = float(fPEl[1])
                            else:
                                fEl[1] = float(fEl[1])*float(fPEl[1])
                            fEls[idx][1] = fEl[1]
                        fList.extend(fEls)
                else:
                    fList.extend([[item[0],int(item[1])] if item[1] != '' else [item[0],1] for item in reFEls.findall(formula)])
                    
                fD = {}
                for k,v in fList:
                    if k in fD:
                        fD[k] = fD[k]+v
                    else:
                        fD[k]=v
                    
                    formula = ''.join(k+str(fD[k]) if fD[k] != 1 else k for k in sorted(fD))
                    
                db[specKey].update({'formula':formula})
    return db
    
if __name__ == "__main__":
    """ 
    Import HbpIS.f90 into a dictionary (not compatible with a mpdb
    dictionary) and then insert key/vals from this dictionary into a
    mpdb compatible dictionary.
    """
    cfg = {}
    cfg['hbpSIf90FileName'] = 'HbpSI.f90'
    cfg['hbpSIf90FQPN'] = os.path.join(mpdb.pref['dataDir'],
                                       'rSander',
                                       cfg['hbpSIf90FileName'])
    rSdb = importHbpSIf90(cfg) # not compatible with mpdb

    """
    buildDT = strftime("%Y%m%d-%H%M")
    db = dict()
    dbmd = dict()
    rSander2mpdb(rSdb,db,dbmd)
    from mpdb.utils import dbSave
    dbSave(db,dbmd,fileName='rSander4mpdb-' + buildDT + '.pickle')
    """
