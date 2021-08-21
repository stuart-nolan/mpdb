"""
mpdb: python Material Property Data Base
slopImport.py: import thermophysical data from the slopYY.dat file 
               (i.e. Gef, Hef, Sef, heat capacity data, etc.)

WWW: https://gitlab.com/ENKI-portal/geopig/-/tree/master/slop
    
Revision Date: 2021.08.15

SPDX-License-Identifier: BSD-2-Clause
Copyright (c) 2021 Stuart Nolan. All rights reserved.
"""
import os
import re
import sys
import pickle
import json
from time import strftime
from pathlib import Path
import mpdb

def iHCP(sdb,name,vals,phaseID,inSec5):
    """
    import Heat Capacity Parameters in sec 1 to 5 to sdb[name]
    """
    values = [float(val) for val in vals]
    sPID = str(phaseID)
    
    if len(values) == 7:
        hCPL = [val for idx,val in enumerate(values) if idx <=2]
        sdb[name].update({'hCP_'+sPID:hCPL})
        sdb[name].update({'TMax_'+sPID:values[3]})
        sdb[name].update({'deltaH_'+sPID:values[4]})
        sdb[name].update({'deltaV_'+sPID:values[5]})
        sdb[name].update({'dPdT_'+sPID:values[6]})
    elif phaseID > 1:
        sdb[name].update({'hCP_'+sPID:values})
    elif inSec5:
        sdb[name].update({'hCP':values})
    else:
        sdb[name].update({'hCP':values})

def iLine4(sdb,name,vals,inSec6):
    """
    import Line 4 of sections 1-5 to db[name]
    """
    values = [float(val) for val in vals]
    sdb[name].update({'deltaG':values[0]})
    sdb[name].update({'deltaH':values[1]})
    sdb[name].update({'entropy':values[2]})
    if not inSec6:
        sdb[name].update({'volume':values[3]})

def importSlop(cfg,sdb,sdbmd):
    """    
    Parameters:
      cfg, configuration dictionary

      cfg['slopFileName'], file name of slopYY.dat file

      cfg['slopFQPN'], fully qualified path name of slopYY.dat file

      sdb, dictionary of dictionaries with the name of components 
           used as top level keys.  Lower level dictionary items 
           are data elements in the slopYY.dat file for each 
           compoennt.
            
      sdbmd, dictionary of dictionaries with meta data for 
             sdb
    """
        
    reDataLine = re.compile(r'^\s.*')
    reComment = re.compile(r'^\*+')
    reFileRev = re.compile(r'^\*+ last updated: (.*?)\s+GEOPIG.*')
    reLine1or2 = re.compile(r'\s*(\S*?)\s+(\S*)\s*')
    reLine3 = re.compile(r'\s*[Rr][Ee][Ff]:(\S*?)\s+(\S*)\s*')
    reValues = re.compile(r'[-+]?\d*\.\d+|[-+]?\d+')
    
    sec1 = "minerals that do not undergo phase transitions"
    inSec1 = False
    sec2 = "minerals that undergo one phase transition"
    inSec2to4 = False
    inSec2 = inSec2to4
    sec3 = "minerals that undergo two phase transitions"
    inSec3 = inSec2to4
    sec4 = "minerals that undergo three phase transitions"
    inSec4 = inSec2to4
    sec5 = "gases"
    inSec5 = False
    sec6 = "aqueous species"
    inSec6 = False
    secLineNum = 0
    lineNum = 0
        
    with open(cfg['slopFQPN'],'r') as dF:
        print("Processing: \"%s\"" % cfg['slopFQPN'])
        cfg['dataSourceDateTime'] = strftime("%Y.%m.%d-%H:%M:%S")
        for line in dF:
            lineNum+=1
            #print('Line Num: %s, %s, %s' % (lineNum,inSec1,inSec2to4))
            if reComment.match(line):
                if reFileRev.match(line):
                    sdbmd['last_updated'] = reFileRev.match(line).groups()[0] 
                continue
            elif reDataLine.match(line) is not None:
                if line.strip() == sec1:
                    inSec1 = True
                    continue
                elif line.strip() == sec2:
                    inSec1 = False
                    inSec2 = True
                    inSec2to4 = True
                    continue
                elif line.strip() == sec3:
                    inSec2 = False
                    inSec3 = True
                    continue
                elif line.strip() == sec4:
                    inSec3 = False
                    inSec4 = True
                    continue
                elif line.strip() == sec5:
                    inSec2to4 = False
                    inSec4 = inSec2to4
                    inSec5 = True
                    continue
                elif line.strip() == sec6:
                    inSec5 = False
                    inSec6 = True
                    continue
                elif line.strip() == '':
                    return
                else:
                    pass
    
            if inSec1 or inSec5:
                secLineNum += 1 #secLineNum should be initialized at 0
                if secLineNum == 1:
                    name = reLine1or2.match(line).groups()[0]
                    sdb.update({name:{},})
                    sdb[name].update({'name':name,
                                     'structural_chemical_formula':
                                     reLine1or2.match(line).groups()[1]})
                elif secLineNum == 2:
                    sdb[name].update({'abbreviation':
                                     reLine1or2.match(line).groups()[0],
                                     'elemental_chemical_formula':
                                     reLine1or2.match(line).groups()[1]})
                elif secLineNum == 3:
                    sdb[name].update({'ref':
                                     reLine3.match(line).groups()[0].split(','),
                                     'date_last_revisited':
                                     reLine3.match(line).groups()[1]})
                elif secLineNum == 4:
                    values = reValues.findall(line)                        
                    if inSec1:
                        sdb[name].update({'state':'solid'})
                    elif inSec5:
                        sdb[name].update({'state':'vapor'})      
                    iLine4(sdb,name,values,inSec6)
                elif secLineNum == 5:
                    values = reValues.findall(line)                        
                    iHCP(sdb,name,values,0,inSec5)
                elif secLineNum == 6:
                    values = reValues.findall(line)
                    if inSec1:
                        sdb[name].update({'TMax':float(values[0])})
                    elif inSec5:
                        sdb[name].update({'TMax':float(values[0])})
                    secLineNum = 0
                    
            elif inSec2to4:
                secLineNum += 1 #secLineNum should be initialized at 0
                    
                if inSec2:
                    lastSecLineNum = 7
                elif inSec3:
                    lastSecLineNum = 8
                elif inSec4:
                    lastSecLineNum = 9
                
                if secLineNum == 1:
                    name = reLine1or2.match(line).groups()[0]
                    sdb.update({name:{}})
                    sdb[name].update({'name':name,
                                     'structural_chemical_formula':
                                     reLine1or2.match(line).groups()[1]})

                elif secLineNum == 2:
                    sdb[name].update({'abbreviation':
                                     reLine1or2.match(line).groups()[0],
                                     'elemental_chemical_formula':
                                     reLine1or2.match(line).groups()[1]})
                elif secLineNum == 3:
                    sdb[name].update({'ref':
                                     reLine3.match(line).groups()[0].split(','),
                                     'date_last_revisited':
                                     reLine3.match(line).groups()[1]})
                elif secLineNum == 4:
                    values = reValues.findall(line)
                    sdb[name].update({'state':'solid'})
                    iLine4(sdb,name,values,inSec6)
                    phaseID = 1 # used to id phase for next set of lines
                elif secLineNum >= 5 and secLineNum < lastSecLineNum:
                    values = reValues.findall(line)
                    iHCP(sdb,name,values,phaseID,inSec5)
                    phaseID += 1
                elif secLineNum == lastSecLineNum:
                    values = reValues.findall(line)
                    sPID = str(phaseID-1)
                    sdb[name].update({'TMax_'+sPID:float(values[0])})
                    secLineNum = 0
                    phaseID = 1
    
            elif inSec6:
                secLineNum += 1 #secLineNum should be initialized at 0
                if secLineNum == 1:
                    name = reLine1or2.match(line).groups()[0]
                    sdb.update({name:{}})
                    sdb[name].update({'name':name,
                                     'structural_chemical_formula':
                                     reLine1or2.match(line).groups()[1]})
                elif secLineNum == 2:
                    sdb[name].update({'abbreviation':
                                     reLine1or2.match(line).groups()[0],
                                     'elemental_chemical_formula':
                                     reLine1or2.match(line).groups()[1]})
                elif secLineNum == 3:
                    sdb[name].update({'ref':
                                     reLine3.match(line).groups()[0].split(','),
                                     'date_last_revisited':
                                     reLine3.match(line).groups()[1]})
                elif secLineNum == 4:
                    values = reValues.findall(line)
                    sdb[name].update({'state':'aqueous'})
                    iLine4(sdb,name,values,inSec6)
                elif secLineNum == 5:
                    values = reValues.findall(line)
                    values = [float(val) for val in values]
                    sdb[name].update({'a1':values[0]})
                    sdb[name].update({'a2':values[1]})
                    sdb[name].update({'a3':values[2]})
                    sdb[name].update({'a4':values[3]})
                elif secLineNum == 6:
                    values = reValues.findall(line)
                    values = [float(val) for val in values]
                    sdb[name].update({'c1':values[0]})
                    sdb[name].update({'c2':values[1]})
                    sdb[name].update({'omega':values[2]})
                    sdb[name].update({'charge':values[3]})                    
                    inSec2to4 = False
                    secLineNum = 0
                    
def saveSlopData(slopDict, buildDT, slopSuffix='MetaData'):
    jsonFileName = 'slop'+ slopSuffix + '-' + buildDT + '.json'
    jsonFQPN = os.path.join(mpdb.pref['dataDir'],jsonFileName)
    with open(jsonFileName,'w') as jsonFile:
        json.dump(slopDict, jsonFile)
        print("Saved slop data to file: \"%s\"" % jsonFQPN)
    
def slopdb2mpdb(sdb,sdbmd,db,dbmd):
    """
    post processing actions to create a "mpdb" compatible pickle 
    file from a slop pickle file.
    """
    from mpdb.periodicTable import fW
    from mpdb.units import uc

    reFormula = re.compile(r'(.*?)\(([0-9]*)\)')
    ref_T = sdbmd['ref_T']
    T_unit = sdbmd['T_unit']
    # unit conversions for cal to J and bar to Pa
    oEU = "cal"
    nEU = "J"
    oPU = "bar"
    nPU = "Pa"
    oVU = "cm^3"
    nVU = "m^3"
    oEU2nEU = uc(1,oEU,nEU)
    oPU2nPU = uc(1,oPU,nPU)
    oVU2nVU = uc(1,oVU,nVU)
    # convert ref P in sdbmd
    ref_P = sdbmd['ref_P']*oPU2nPU
    ref_P_unit = sdbmd['ref_P_unit'].replace(oPU,nPU)
    
    # create converstion Factors (cFs) from unit conversions and scaling
    # factors
    hCP_cFs = []
    hCP_units = {}
    for key in ['a_','b_','c_']:
        hCP_cFs.append(oEU2nEU/sdbmd[key + 'scale'])
        hCP_units.update({key + 'unit':sdbmd[key + 'unit'].replace(oEU,nEU)})
    hC_unit = sdbmd['a_unit'].replace(oEU,nEU)

    xi_cF = {}
    xi_unit = {}
    for key in ['a1', 'a2', 'a3', 'a4', 'c1', 'c2']:
        if oEU in sdbmd[key + '_unit']:
            eUC = oEU2nEU
        else:
            eUC = 1
        if oPU in sdbmd[key + '_unit']:
            pUC = oPU2nPU
        else:
            pUC = 1
            
        xi_cF[key] = eUC/pUC/sdbmd[key + '_scale']
        xi_unit.update({key:sdbmd[key + '_unit'].replace(oEU,nEU)})
        xi_unit.update({key:xi_unit[key].replace(oPU,nPU)})

    # flags for doing unit conversions not already done above
    for key in sdb.keys():
        if key not in db.keys():
            db.update({key:{}})
            db[key].update({'name':key})
            
        dbmd.update({key:{}})
        dbmd[key].update({'slopFile':sdbmd['slopFile']})
        dbmd[key].update({'rev':sdbmd['last_updated']})
        dbmd[key]['ref'] = sdb[key]['ref']
        dbmd[key]['date_last_revisited'] = sdb[key]['date_last_revisited']
        #
        # BEGIN formula and fW update to slop values
        #
        formList = reFormula.findall(sdb[key]['elemental_chemical_formula'])
        formula = ''
        fWeight = 0
        for item in formList:
            if item[0] == '+' or item[0] == '-':
                continue
            elif item[1] != '1':
                formula += item[0]+item[1]
            else:
                formula += item[0]
            
            fWeight += fW(item[0])*float(item[1])
            
        db[key].update({'formula':formula})
        db[key].update({'fW':fWeight})
        db[key].update({'state':sdb[key]['state']})
        #
        # END add formula and fW key & val
        #
        #
        # BEGIN rescale and convert energy, pressure & volume units
        #
        kks = sdb[key].keys()
        if 'charge' in kks:
            db[key].update({'charge':sdb[key]['charge']})
            
        for item in ['hCP','hCP_1','hCP_2','hCP_3','hCP_4']:
            if item in kks:
                db[key].update({item:[val*hCP_cFs[idx] for idx,val in
                                          enumerate(sdb[key][item])]})
                dbmd[key].update({'T_unit':T_unit})
                dbmd[key].update({'hC_unit':hC_unit,
                                        'hCP_description':
                                        'heat Capacity Parameters: hC = hCP[0] + hCP[1]*T + hCP[2]/T^2, where T is in K'})
        for item in ['deltaG','deltaH','entropy','omega']:
            if item in kks:
                if item == 'deltaG':
                    newItem ='eFG298K'
                    val = sdb[key][item]
                    if val != 999999.0:
                        db[key].update({newItem:oEU2nEU*val/1000})
                        desc = 'Gibbs Free Energy of Formation at %s %s and %s %s' % (ref_T,T_unit,ref_P,ref_P_unit)
                        uKey = newItem + '_unit'
                        dbmd[key].update({uKey:"k{0}".format(sdbmd[item + '_unit'].replace(oEU,nEU)),
                                            newItem + '_description':desc})
                if item == 'deltaH':
                    newItem ='eFH298K'
                    val = sdb[key][item]
                    if val != 999999.0:
                        db[key].update({newItem:oEU2nEU*val/1000})
                        desc = 'Enthalpy of Formation at %s %s and %s %s' % (ref_T,T_unit,ref_P,ref_P_unit)
                        uKey = newItem + '_unit'
                        dbmd[key].update({uKey:"k{0}".format(sdbmd[item + '_unit'].replace(oEU,nEU)),
                                        newItem + '_description':desc})
                if item == 'entropy':
                    newItem = 'eFS298K'
                    val = sdb[key][item]
                    if val != 999999.0:
                        db[key].update({newItem:oEU2nEU*val})
                        desc = 'Entropy of Formation at %s %s and %s %s' % (ref_T,T_unit, ref_P,ref_P_unit)
                        dbmd[key].update({newItem+'_unit':sdbmd[item + '_unit'].replace(oEU,nEU),
                                            newItem+'_description':desc})
                if item == 'omega':
                    db[key].update({item:oEU2nEU*sdb[key][item]})
                    desc = 'Reference Born coefficient for the \"HKF\" EOS'
                    dbmd[key].update({item+'_unit':sdbmd[item + '_unit'].replace(oEU,nEU),
                                        item+'_description':desc})
                    # rescale omega
                    db[key][item] = db[key][item]/sdbmd[item + '_scale']
                    
        for item in ['deltaH_1','deltaH_2','deltaH_3','deltaH_4']:
            if item in kks:
                val = sdb[key][item]
                if val != 999999.0:
                    db[key].update({item:val*oEU2nEU})

                    desc = 'Enthalpy change at x to x+1 phase transition T & P'
                    item = 'deltaH_x'

                    dbmd[key].update({'deltaH_unit':sdbmd['deltaH_unit'].replace(oEU,nEU),
                                        item + '_description':desc})

        for item in ['TMax','TMax_1','TMax_2','TMax_3','TMax_4']:
            if item in kks:
                val = sdb[key][item]
                if val != 999999.0:
                    db[key].update({item:val})
                    if item == 'TMax':
                        desc = 'Max temperature for heat capacity parameters (hCP)'
                    else:
                        desc = 'Max temperature for phase x heat capacity parameters (hCP_x)'

                    dbmd[key].update({'T_unit':sdbmd['T_unit'],
                                        'TMax_description':desc})

        for item in ['dPdT','dPdT_1','dPdT_2','dPdT_3','dPdT_4']:
            if item in kks:
                val = sdb[key][item]
                if val != 999999.0:
                    db[key].update({item:val*oPU2nPU})
                    # add units, descriptions to dbmd[key]
                    if item == 'dPdT':
                        desc = 'dP/dT at %s %s and %s %s' % (ref_T,T_unit,ref_P,ref_P_unit)
                    else:
                        desc = 'dP/dT at x phase transition T & P'
                        item = 'dPdT_x'
                    dbmd[key].update({'dPdT_unit':sdbmd['dPdT_unit'].replace(oPU,nPU),
                                        item+'_description':desc})
                    
        for item in ['volume','deltaV_1','deltaV_2','deltaV_3','deltaV_4']:
            if item in kks:
                val = sdb[key][item]
                if val != 999999.0:
                    db[key].update({item:val*oVU2nVU})

                    if item == 'volume':
                        desc = 'Volume at %s %s and %s %s' % (ref_T,T_unit,ref_P,ref_P_unit)
                    else:
                        desc = 'delta volume at x phase transition T & P'
                        item = 'deltaV_x'

                    dbmd[key].update({'volume_unit':sdbmd['volume_unit'].replace(oVU,nVU),
                                        item + '_description':desc})
        
        for item in ['a1', 'a2', 'a3', 'a4', 'c1', 'c2']:
            if item in kks:
                db[key].update({item:sdb[key][item]*xi_cF[item]})
                desc = 'Regression coefficients used by the \"HKF\" EOS'
                dbmd[key].update({item+'_unit':xi_unit[item],
                                    'ai_ci_description':desc})
        #
        # END rescale and convert values to SI units
        #
    
if __name__ == "__main__":
    """ 
    Import the data and meta data from the slopYY.dat file into the
    "sdb" (data) and the "sdbmd" (meta data) dictionaries.  Convert
    and insert items from the sdb and sdbmd dictionaires into their
    respective "mpdb" compatible dictionries.
        
    slop16.dat, rev. 29.May.19:
    line ~1184 (component, TITANITE) add "ref:" in front of "7,9,19"
    line ~2044 (component: O2,g) change o(2) to O(2)
    change lines 3594 and 3595 from:
        H+                  h(+)
        H+                  h(1)+(1)        
    to:
        H+                  H(+)
        H+                  H(1)+(1)
    change line 9628, value "14 .9410" to "14.9410"
    """
    path = os.path.join(str(Path.home()),
                        'local',
                        'src',
                        'gitlab',
                        'geopig',
                        'slop')
    cfg = {}
    cfg['slopFileName'] = "slop16.dat"
    cfg['slopFQPN'] = os.path.join(path, cfg['slopFileName'])
    sdb = dict()
    sdbmd = dict()
    # import slop db meta data dictionary from json file
    restoreDT = '20210814-1950'
    with open('slopMetaData-' + restoreDT + '.json') as jsonFile:
        sdbmd = json.load(jsonFile)
    
    buildDT = strftime("%Y%m%d-%H%M")
    # uncomment to save the (meta) data import to a json file
    saveSlopData(sdbmd, buildDT, slopSuffix='MetaData')

    importSlop(cfg,sdb,sdbmd)
    saveSlopData(sdb, buildDT, slopSuffix='Data')

    #
    # Post processing to convert and insert sdb and sdbmd into "mpdb"
    # compatible dictionaries
    #
    db = dict()
    dbmd = dict()
    slopdb2mpdb(sdb,sdbmd,db,dbmd)
    from mpdb.utils import dbSave
    dbSave(db,dbmd,fileName='slop4mpdb-' + buildDT + '.pickle')
