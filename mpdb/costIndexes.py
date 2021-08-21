"""mpdb: Material Property Data Base as python package
costIndexes.py: fetch annual average cost indexes from 
                data/costIndexes/costIndex.csv

Revision Date: 2021.08.18

SPDX-License-Identifier: BSD-2-Clause
Copyright (c) 2021 Stuart Nolan. All rights reserved.
"""
import os
import pathlib
import csv
import mpdb

def cepci(year, ciKey="ce plant cost index", idxKey="year"):
    """Wrapper function for costIndex to fetch the annual average
    chemical engineering plant cost index for the requested year.

    pyUNO, libreOffice MasterScriptProviderFactory notes:

    When python scripts are invoked from libreOffice basic, the python
    __file__ variable is overwrtitten by pyUNO, sys.argv[0] is empty,
    sys.path[0] is "/usr/lib/libreoffice/program", os.getcwd() and
    pathlib.Path().resolve() return the user's home directory.

    To install, link to this file in the libreoffice script directory: 

    ubuntu: ~/.config/libreoffice/4/user/Scripts/python 

    then create a pyUno compatible libreOffice script to call cepci
    below as usual and set mpdb.pref['costIndexFile'] to the fully
    qualified path name of costIndexes.csv.

    Parameters: 
    - year, YYYY, i.e. "1997"
    - ciKey, cost index key in the data dictionary
    - idxKey, index key in the data dictionary used for ciKey values

    Returns:
    - annual average cost index for year as a float or
      None if no data is found.

    """
    cif = mpdb.pref['costIndexFile']
    ci = costIndex(cif)
    return ci.getIndex("%s" % year, ciKey=ciKey, idxKey=idxKey)

def cepci_keys():
    cif = mpdb.pref['costIndexFile']
    ci = costIndex(cif)
    keys = ""
    for key in ci.data.keys():
        keys =  keys + key + ','
    return keys[:-1]
    
class costIndex():
    def __init__(self,dataFile):
        self.dataFQPN = dataFile
        self.data = self.importCSV()
    
    def getIndex(self, year, ciKey='ce plant cost index', idxKey='year'):
        """
        fetch annual average cost index for year as a float

        Parameters: 
        - year, YYYY, i.e. "1997"
        - ciKey, cost index key in the data dictionary
        - idxKey, index key in the data dictionary used for ciKey values

        Returns:
        -  if it exists, annual average cost index for year as a float.
           None otherwise.
        """
        if idxKey in self.data.keys():
            if year not in self.data[idxKey]:
                return None
            ci = self.data[ciKey][self.data[idxKey].index(year)]
            if ci != '':
                return float(ci)
            else:
                return None
        else:
            return self.data[__file__]

    def importCSV(self):
        """
        python open function requires full path to file if the file 
        is not in the script working directory
        """    
        try:
            csvFile = open(self.dataFQPN, 'r')
            reader = csv.reader(csvFile,dialect='excel')
        except:
            errMsg = "File: \"%s\" not found" % (self.dataFQPN)
            return { __file__ : errMsg }

        rows = []
        header = next(reader,None)
        data = dict()
        for h in header:
            data[h] = []
        for row in reader:        
            for h, v in zip(header, row):
                data[h].append(v)
            
        csvFile.close()
        return data

if __name__ == "__main__":
    # Import
    mpdb.pref['costIndexFile'] = os.path.join(mpdb.pref['dataDir'],
                                              'costIndexes',
                                              'costIndexes.csv')
    mpdb.savePref(pref=mpdb.pref)
    cif = mpdb.pref['costIndexFile']
    ci = costIndex(cif)

    # Test
    print(cepci_keys())
    year = '2020'
    print("CEPCI in %s: %s" % (year,cepci(year)))
    ciKey='ppi industrial chemicals'
    print("\"%s\" in %s: %s" % (ciKey,year,ci.getIndex(year, ciKey=ciKey)))
