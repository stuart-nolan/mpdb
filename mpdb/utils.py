"""
mpdb: Material Property Data Base as python module
utils.py: mpdb utilities

Revision Date: 2021.08.15

SPDX-License-Identifier: BSD-2-Clause
Copyright (c) 2021 Stuart Nolan. All rights reserved.

1. The data:

   YCPH data are csv files that are imported into a class instance
   with a "mpdb" compatilbe dictionary.  "mpdb" attributes are stored
   in a separate "mpdbmd" (mpdbMetaData) dictionary contained in the
   same class instance.

   rSander data are static f90 files that are imported into python
   dictionaries (incompatible with "mpdb/mpdbmd" dictionaires) and then
   converted into dictionaries compatible with "mpdb/mpdbmd" dictionaries.
          
   slop data source is slopYY.dat which is imported into python
   dictionaries and then converted into dictionaries compatible with
   "mpdb/mpdbmd" dictionaries.

   - edit data in the original data sources
   - recreate dictionary & pickle after edditing data
   - limited versioning for (potentially edited) original data sources
   - share excerpts from YAWS and CEPCI?
   - export mpdb dictionaies to csv files and/or sqlite db

2. need to limit project scope and turn it into a python (setup.py) package
  - is cython really needed (esp. if used from spreadsheets)?
  - one top level function (__init__.py?) to load units, dbLoad, 
    dbFind, periodicTable
  - function to find dependencies like gnu units exe
  - directory structure: data, eq, gui, utils (for inital data source 
    imports), lo (libreoffice), excel, 

3. (python) GUI
   - Tkinter too limited?
     - ttkFilterTree loads large (5000 components) data sets too slowly
     - only load ~500 common componets but still be able to search all?
   - Qt, WxPython (GTK based, not actively developed?) have been tried.  
   - Kivy?
   - means (GUI) to update and edit data in mpdb dictionaires

4. integration with spreadsheets
   - excel, excel online
   - libreoffice
   - google sheets

5. Calculation functions
   - generalized to call calculation functions (e.g. CoolProp, YAWS 
     equations, etc.) based on the property (i.e. density, formula weight, 
     etc.) and not the data source
   - create VLE database for NRTL
   - NRTL to cython
   - create gas interactithermoFuncson database for pr,srk
   - simple multicomponent streams & material balance (mix, recycle, etc)
   - basic reaction energy handling
   - basic phase transistion handling
   - basic multicomponent from pure component calculations

6. units, bounds, 
   - convert qunatities to SI units
   - equations and values need bounds
     - pint, units, physics, scimath,

"""

import os
from pathlib import Path
import re
import itertools as IT
from copy import deepcopy
from time import strftime
import pickle
from tabulate import tabulate
import mpdb

def dbFind(sVal,sKeys=['name','formula'],rKeys=['formula','CAS'],db={}):
    """
    dbFind(sVal,sKeys=['name','formula'],rKeys=['formula','CAS'],db={})
    
    search in db[key] for all key in db.keys()
    
    Paramters:
        sVal, search string (a python regex, e.g. .*Al.*)
        sKeys, list of sKey's to search in db[key] i.e. db[key][sKey] (for all sKey in sKeys)
        rKeys, list of rKey's in db[key] to return i.e. db[key][rKey] (for all rKey in rKeys)
        db, db dictionary instance
        
    Returns:
        list of lists in the form [[index, key, db[key][rKey], ...] for each
        rKey in rKeys and for each sKey in sKeys for which db[item][sKey] 
        matches 'sVal'

    Notes:
        normally key in db[key] are material names in aqueous or yaws db's
    """
    idx = IT.count(1)
    res = []
    
    regex = r'%s' % sVal
    try:
        fKey = re.compile(regex)
    except:
        return print("val: %s is not a valid regex" % sVal)

    for k in db.keys():
        for sKey in sKeys:
            if sKey in db[k].keys():
                if fKey.match(db[k][sKey]):
                    resEl = [next(idx),k]
                    resEl.extend([db[k][rKey] for rKey in rKeys if rKey in db[k].keys()])
                    res.append(resEl)
                    break #prevent duplicate entries
    return res

def dbKeys(db):
    """
    Fetch name, CAS, and formula from db and return them as a list of tuples

    Parameters:
        db, nested dictionary like {'<name>':{'CAS':'', 'formula':'', ...}}
    
    Returns:
        nCF, list of tuples like [('name0','CAS0','formula0'),
                                  ('name1','CAS1','formula1'),
                                  ...]
    Notes:
        extract names, CAS, and formulas to separate tuples via zip, e.g.:
            dbNames, dbCAS, dbFormula = zip(*dbKeys(db))
        
    """
    nCF = []
    for k in db.keys():
        if 'CAS' in db[k].keys():
            CAS = db[k]['CAS']
        else:
            CAS = ""
        if 'formula' in db[k].keys():
            formula = db[k]['formula']
        else:
            formula = ""
        nCF.append((k,CAS,formula))
        
    return nCF

def dbLoad(fileName='mpdb-YYYYMMDD-HHMM.pickle',
           interactive=False,
           dataDir=mpdb.pref['dataDir']):
    """
    dbLoad(fileName='mpdb-YYYYMMDD-HHMM.pickle',
           interactive=False,
           dataDir=mpdb.pref['dataDir'])
    
    Load a pickle file saved by dbSave
    
    Parameters:
        fileName, pickle file to load.  Must include the path to the file 
                  in fileName.  dataDir is ignored.
        interactive, True to use text base file selection in dataDir
        dataDir, root directory containing "*.pickle" files.  Only valid
                 if interactive=True.  Only files ending in ".pickle" are
                 listed.
    
    Returns:
        (db, dbmd)
    """
    if interactive:
        #read all pickle files
        
        fileNames = []
        for (root, dirs, files) in os.walk(Path(dataDir)):
            for file in files:
                if file.endswith(".pickle"):
                    fileNames.append(os.path.join(root,file))

        nFiles = len(fileNames)
        header = ['Index','File Name']
        inStr = "Select an index (%s to %s) or e(x)it: " % (1,nFiles)
        choice = ''
        while choice != 'x':
            table = []
            for idx in range(nFiles):
                table.append([idx+1])
                table[idx].append(fileNames[idx])
            print('\n'+tabulate(table,headers=header,numalign='center',
                                stralign='left')+'\n')
            choice = input(inStr)
            invalStr = "Invalid input: %s\n0 <= index <= %s or \'x\' to exit\n" % (choice,nFiles)
            if re.match(r'^[0-9]+$',choice):
                index = int(choice)
                if index > nFiles: #check for valid index
                    inStr = invalStr + inStr
                else:
                    fileName = fileNames[index-1]
                    fileName = os.path.join(fileName)
                    print(fileName)
                    break 
            else:
                inStr = invalStr + inStr

    if fileName == 'mpdb-YYYYMMDD-HHMM.pickle':
        return ({},{})
        
    with open(fileName, 'rb') as f:
        return pickle.load(f)

def dbMerge_byCAS(db,dbmd,idb,idbmd):
    """
    merge idb and idbmd into db and dbmd by CAS.  If
    db['<name>']['CAS'] or if idb['<inputName>']['CAS'] do not exist,
    then items from idb/idbmd are not merged into db.
    
    Parameters:
        db, dictionary like {'<name>':{'CAS':'', ...}}
        dbmd, dictionary like {'<name>':{...}}
        idb, dictionary like {'<inputName>':{'CAS':'', ...}}
        idbmd, dictionary like {'<inputName>':{...}}
    
    Notes:
        - '<name>' key must exist and be the same in db and dbmd
        - '<inputName>' key key must exist and be the same in idb and 
          idbmd
        - '<name>' and '<inputName>' need not match

    """
    dbNames, dbCAS, dbFormula = zip(*dbKeys(db))
    
    for key,val in idb.items():
        if 'CAS' in val.keys():            
            if val['CAS'] in dbCAS and val['CAS'] != "":
                dbName = dbNames[dbCAS.index(val['CAS'])]
                for k,v in val.items():
                    if k not in db[dbName].keys():
                        db[dbName].update({k:v})
                for k,v in idbmd[key].items():
                    if k not in dbmd[dbName].keys():
                        dbmd[dbName].update({k:v})

def dbMerge_byName(db,dbmd,idb,idbmd,preserve=False):
    """
    merge idb and idbmd into db and ibdmd by '<name>'.  If
    idb['<name>'] exits in db['<name>'] then the val of idb['<name>']
    is updated in db['<name>'].  Otherwise, idb['<name>'] is updated
    in db.
    
    Parameters:
        db, dictionary to be updated
        dbmd, dictionary of meta to be updated
        idb, input dictionary/data to be merged into db
        idbmd, input dictionary/meta data to merged into dbmd
        preserve - boolean, if true idb[key] is not added to db[key] if key
                   exits in both db and idb

    """
    for k in idb.keys():
        if k in db.keys():
            if preserve == True:
                pass
            else:
                db[k].update(idb[k])
                if k in idbmd.keys():
                    for key,val in idbmd.items():
                        if dbmd[k][key] in idbmd[k].keys():
                            dbmd[k][key].update(idbmd[k])
        else:
            db.update({k:idb[k]})
            dbmd.update({k:idbmd[k]})

    return (db,dbmd)

def dbSave(db,
           dbmd,
           fileName='mpdb-' + strftime("%Y%m%d-%H%M") + '.pickle'):
    """
    dbSave(db,
           dbmd,
           fileName=mpdb-' + strftime("%Y%m%d-%H%M") + '.pickle'):

        save a pickle file of [db, dbmd]
    
    Parameters:
        db, dictionary of material properties
        dbmd, dictionary of db properties  
        fileName, pickle file FQPN
          - fileName is autmatically added to dbmd['pickleFQPN']
    Returns:
        fileName: pickle File Fully Qualified Path Name
    """
   
    if isinstance(dbmd,dict):
        dbmd.update({'pickleFQPN':fileName})
        
    with open(fileName, 'wb') as f:
        pickle.dump([db,dbmd], f)
    
    return fileName
    
def dbSubset(keys, db):
    """
    dbSubset(keys, db)
    
    create a deepcopy of db for each key in keys
    
    Parameters:
        keys, (partial) list of keys in db
        db, material property database dictionary instance
    
    Returns:
        dbSubset, deepcopy of db for each key in keys
    """
    dbSubset = deepcopy({key:val for key,val in db.items() if key in keys})
    return dbSubset

def rSanderSelect(dbItem,index=0,interactive=False):
    """
    rSanderSelect(dbItem,index=0,interactive=False)
        
    select which rSander henry data to use in dbItem
    
    Parameters:
        dbItem, db[key] dictionary object with keys = ['hbpSIP','hbpSIPL',
                'hbpSI_index']
        index, positive integer index for list item in hbpSIPL to move into 
               hbpSIP.  Use interactive=True to display choices and ask for 
               user input for an index.
        interactive, True to display choices and ask user input for an index, 
                     False to make the change silently
    Returns:
        Nothing on success (dbItem is succussfully changed) or error messages
        if there is an issue
    """
    keys = ['hbpSIP','hbpSIPL','hbpSI_index']
    for key in keys: #test if dbItem has valid dictionary keys
        if key not in dbItem.keys():
            return print("Henry data (%s) not found in dbItem[%s]\n" % 
                         (key,dbItem['name']))
    nHbpSIPL =len(dbItem['hbpSIPL'])
    if not interactive:
        invalIndex = "Invalid index: %s\n0 <= index <= %s\n" % (index,nHbpSIPL-1)
        if re.match(r'^[0-9]+$',str(index)):  #make sure index is positive integer
            if index > nHbpSIPL-1: #check for valid index
                return print(invalIndex)
            dbItem['hbpSI_index'] = index
            dbItem['hbpSIP'] = [float(dbItem['hbpSIPL'][index][0]),
                                float(dbItem['hbpSIPL'][index][1])]
        else:
            return print(invalIndex)
    else:
        header = ['Index','Ho /mol/kg/Pa','dln(H)/d(1/T) /K','Code','Ref.']
        inStr = "Select an index (%s to %s) or e(x)it: " % (0,nHbpSIPL-1)
        choice = ''
        while choice != 'x':
            table = []
            for idx in range(nHbpSIPL):
                table.append([idx])
                table[idx].extend(dbItem['hbpSIPL'][idx])
            print('\n'+tabulate(table,headers=header,numalign='center',
                                stralign='center')+'\n')
            choice = input(inStr)
            inStr = "Select an index (%s to %s) or e(x)it: " % (0,nHbpSIPL-1)
            invalStr = "Invalid input: %s\n0 <= index <= %s or \'x\' to exit\n" % (choice,nHbpSIPL-1)
            if re.match(r'^[0-9]+$',choice):
                index = int(choice)
                if index > nHbpSIPL-1: #check for valid index
                    inStr = invalStr + inStr
                else:
                    dbItem['hbpSI_index'] = index
                    dbItem['hbpSIP'] = [float(dbItem['hbpSIPL'][index][0]),
                                        float(dbItem['hbpSIPL'][index][1])]
            else:
                inStr = invalStr + inStr

if __name__ == "__main__":
    #%run ycph
    yDataDir = os.path.join(mpdb.pref['dataDir'],'yaws',)
    (db, dbmd) = dbLoad(interactive=True,dataDir=yDataDir)
    
    rSDataDir = os.path.join(mpdb.pref['dataDir'],'rSander',)
    (rSdb, rSdbmd) = dbLoad(interactive=True,dataDir=rSDataDir)

    sDataDir = os.path.join(mpdb.pref['dataDir'],'slop',)
    (sdb, sdbmd) = dbLoad(interactive=True,dataDir=sDataDir)

    dbMerge_byCAS(db,dbmd,rSdb,rSdbmd)
    dbMerge_byName(db,dbmd,sdb,sdbmd)
    
    mpdbFQPN = os.path.join(mpdb.pref['dataDir'],'mpdb-' + 
                            strftime("%Y%m%d-%H%M") + 
                            '.pickle')
    
    #dbSave adds key pickleFQPN to dbmd
    mpdbFileName = dbSave(db,dbmd,fileName=mpdbFQPN)
    print("Created mpdb pickle: %s" % mpdbFileName)
    """
    try:
        if not isinstance(db,dict):
            (db, dbmd) = dbLoad(interactive=True)
    except:
        (db, dbmd) = dbLoad(interactive=True)
    """
