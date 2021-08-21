"""
mpdb: Material Property Data Base as python package
yawsImport.py: Yaws Chemical Property Handbook 

Revision Date: 2021.08.16

SPDX-License-Identifier: BSD-2-Clause
Copyright (c) 2021 Stuart Nolan. All rights reserved.
"""
from mpdb.utils import dbSave, dbLoad
from mpdb.eq.yaws_python import *
import csv
import os
import re
from time import strftime, sleep
from math import log10
import requests
from html.parser import HTMLParser

class pHTML(HTMLParser):
    def __init__(self):
        HTMLParser.__init__(self)
        self.out = []

    def handle_data(self,data):
        self.out.append(data)

def updateCAS():
    """
    find missing CAS numbers and add them in yaws csv files
    
    EDITS MAY BE REQUIRED BEFORE USE
    """
    # import yaws db only...
    yawsPath = os.path.join('yaws')
    (db, dbMetaData) = dbLoad(interactive=True, dataDir=yawsPath)
    noCAS = {key:val['formula'] for (key,val) in db.items() if val['CAS'] == ''}
    URL = "https://toxgate.nlm.nih.gov/cgi-bin/sis/search2"
    PARAMS = {"queryxxx":"", "database":"hsdb", "Stemming":1, "and":1, "second_search":1, "gateway":1, "chemsyn":1}
    p = pHTML()

    CAS = {}
    for key,val in noCAS.items():
        PARAMS["queryxxx"] = key
        print("Looking up %s cas number" % key)
        r = requests.post(url = URL, data = PARAMS)
        p.feed(r.text)
        if 'CAS Registry Number: ' in p.out:
            CAS[key] = p.out[p.out.index('CAS Registry Number: ')+1]
            print("Found %s cas: %s" % (key,CAS[key]))
        p.out = []

        # don't spam the server...
        sleep(1)

    """
    for key,val in CAS.items():
        noCAS.pop(key,None)
    """
    
    fileNames = ['activatedCarbonAdsorption.csv',
                 'criticalProperty.csv',
                 'energyFormationGAS.csv',
                 'energyFormationHUS.csv',
                 'enthalpyCombustion.csv',
                 'enthalpyFusion.csv',
                 'enthalpyVaporization.csv',
                 'henrysLaw.csv',
                 'liquidHeatCapacity.csv',
                 'liquidDensity.csv',
                 'liquidThermalConductivity.csv',
                 'liquidViscosity.csv',
                 'saltWaterSolubility.csv',
                 'solidHeatCapacity.csv',
                 'surfaceTension.csv',
                 'thermalConductivity.csv',
                 'vaporEnthalpyFormation.csv',
                 'vaporGibbsEnergyFormation.csv',
                 'vaporHeatCapacity.csv',
                 'vaporPressure.csv',
                 'vaporThermalConductivity.csv',
                 'vaporViscosity.csv',
                 'waterSolubility.csv']
    #read the data in
    path = 'data'
    for fileName in fileNames:
        try:
            file = os.path.join(path,fileName)
            csvFileHandle = open(file, 'r')
            csvReader = csv.reader(csvFileHandle,dialect='excel')
        except:
            print("File: %s" %file+" not found.  Check path and file name.")
            return        
    
        csvRows = []
        # custom update code
        for row in csvReader:        
            if row[1] in CAS.keys(): #row[1] is the chemical name
                row[2] = CAS[row[1]] #replace CAS number
            csvRows.append(row)
        csvFileHandle.close()
        # rename the original csv file before saving the updated one
        os.rename(file,file+'.bak') 
        
        #output a new csvFile
        csvFileHandle = open(file, 'w')
        csvWriter = csv.writer(csvFileHandle,dialect='excel')
        csvWriter.writerows(csvRows)
        csvFileHandle.close()

class YCPH():
    def __init__(self):
        """
        Pure component chemical substance data from Yaws' Handbook of 
        Thermodynamic and Physical Properties of Chemical Compounds.  
        """
        self.refndSpace = re.compile(r'\s+')
        self.reUnitSep = re.compile(r'(.*?)/(.*)')
        self.reCSV = re.compile(r'.*csv')
        self.csvPath = os.path.join('yaws')
        self.db = {}
        self.dbMetaData = {}
        self.importDateTime = strftime("%Y%m%d-%H%M")
        
    def importYaws(self):
        """
        importYaws() 
            Call the YCPH methods in the specifed sequence.  

        Returns:
            [self.db, self.dbMetaData]
        """
        # initilize db and dbMetaData in criticalProperty
        self.criticalProperty()
        
        # call the following functions after criticalProperty but otheriwse
        # without regard to sequence
        self.activatedCarbonAdsorption()
        self.energyFormationGAS()
        self.energyFormationHUS()
        self.enthalpyFusion()
        self.enthalpyVaporization()
        self.liquidDensity()
        self.liquidHeatCapacity()
        self.liquidThermalConductivity()
        self.liquidViscosity()
        self.saltWaterSolubility()
        self.solidHeatCapacity()
        self.surfaceTension()
        self.thermalConductivity()
        self.vaporEnthalpyFormation()
        self.vaporGibbsEnergyFormation()
        self.vaporHeatCapacity()
        self.vaporPressure()
        self.vaporThermalConductivity()
        self.vaporViscosity()
        self.waterSolubility()

        return[self.db, self.dbMetaData]
        
    def activatedCarbonAdsorption(self):
        """
        see doc strings for criticalProperty
        """
        csvFQPN = os.path.join(self.csvPath,'activatedCarbonAdsorption.csv')
        self.csvFile(csvFQPN)
        header = next(self.csvReader)
        header = [re.sub(self.refndSpace,'',h) for h in header]

        csvDict = {'aCAP':
                       {'colIdx':[3,4,5],
                        'units':'g/100 g',
                        'description':'activated Carbon Adsorption Parameters',
                        'eqnDoc':activatedCarbonAdsorption.__doc__},
                   'aCACmin':
                       {'colIdx':[6],
                        'units':self.reUnitSep.match(header[6]).groups()[1],
                        'description':'activated Carbon Adsorption Concentration Min'},
                   'aCACmax':
                       {'colIdx':[7],
                        'units':self.reUnitSep.match(header[7]).groups()[1],
                        'description':'activated Carbon Adsorption Concentration Max'}}
        #update csvDict for entries common to all top level keys
        for key in csvDict.keys():
            csvDict[key].update({'importDateTime':self.importDateTime,
                                 'dataSource':csvFQPN})

        metaDataKeys = ['description','dataSource','importDateTime',
                        'units','eqnDoc']
        self.csv2db(csvDict,metaDataKeys)
        self.csvFileHandle.close() #close the file openend by csvFile
    
    def criticalProperty(self):
        """
        criticalProperty.csv should be processed first to populate self.db and 
        self.dbMetaData top level dictionary keys with material names.
        """
        #1. open the csv data file
        csvFQPN = os.path.join(self.csvPath,'criticalProperty.csv')
        self.csvFile(csvFQPN)
        #2. get csv header (1st line) and remove white space from header entries
        header = next(self.csvReader)
        header = [re.sub(self.refndSpace,'',h) for h in header]

        #3. define a dictionary to decode the csv file in self.csv2db()
        
        csvDict = {'formula':
                       {'colIdx':[0],
                        'units':None,
                        'description':'elemental formula'
                        },
                   'name':
                       {'colIdx':[1],
                        'units':None,
                        'description':'chemical substance name'
                        },
                    'CAS':
                       {'colIdx':[2],
                        'units':None,
                        'description':'registry number'
                        },
                    'fW':
                       {'colIdx':[3],
                        'units':None,
                        'description':'formula Weight'
                        },
                    'nFP':
                       {'colIdx':[4],
                        'units':self.reUnitSep.match(header[4]).groups()[1],
                        'description':'normal Freezing Point'
                        },
                    'nBP':
                       {'colIdx':[5],
                        'units':self.reUnitSep.match(header[5]).groups()[1],
                        'description':'normal Boiling Point'
                        },
                    'critT':
                       {'colIdx':[6],
                        'units':self.reUnitSep.match(header[6]).groups()[1],
                        'description':'critical Temperature'
                        },
                    'critP':
                       {'colIdx':[7],
                        'units':self.reUnitSep.match(header[7]).groups()[1],
                        'description':'critical Pressure'
                        },
                    'critV':
                       {'colIdx':[8],
                        'units':self.reUnitSep.match(header[8]).groups()[1],
                        'description':'critical Volume'
                        },
                    'critD':
                       {'colIdx':[9],
                        'units':self.reUnitSep.match(header[9]).groups()[1],
                        'description':'critical Density'
                        },
                    'critC':
                       {'colIdx':[10],
                        'units':None,
                        'description':'critical Compressibility'
                        },
                    'acentric':
                       {'colIdx':[11],
                        'units':None,
                        'description':'acentric factor'
                        }}
        #update csvDict for entries common to all top level keys
        for key in csvDict.keys():
            csvDict[key].update({'importDateTime':self.importDateTime,
                                 'dataSource':csvFQPN})
                                  
        metaDataKeys = ['description','dataSource','importDateTime',
                        'units','eqnDoc']
        self.csv2db(csvDict,metaDataKeys)
        self.csvFileHandle.close() #close the file openend by csvFile

    def csvFile(self,csvFQPN):
        """
        utility function to open the requested csv file in the csv dir and 
        create a csv reader for this file
        
        the calling entity should call self.csvFileHandle.close() when done
        with the requested csv file
        
        Parameters:
          csvFQPN, csv Fully Qualified Path Name as string
        """
        try:
            self.csvFileHandle = open(csvFQPN, 'r')
        except:
            print("File: %s not found.  Check path and file name." % csvFQPN)
        
        self.csvReader = csv.reader(self.csvFileHandle,dialect='excel')

    def csv2db(self,csvDict,metaDataKeys):
        """
        Parse rows from a csv file containing yaws chemical property handbook
        data.  For each csv file processed in which the first line contains 
        column header information, next(self.csvReader) should be 
        called once before calling this function.
        
        Parameters
            csvDict: dictionary of the form:
                     {'<dbProp>':{'colIdx':[integer(s)], ...}, ...}  
                     - '<dbProp>' can be any string value used to identify a 
                       property (column) from the csv file input into db - e.g.
                       'formula', 'CAS', 'fW', etc.
                     - Every '<dbProp>' must have a dictionary value with 
                       'colIdx' as a key.  The value for 'colIdx' must be a 
                       list of integers used to identify the column(s) to 
                       extract data from the csv file per '<dbProp>'
                     - If the 'colIdx' value is a list with one integer, the
                       data from this column csv will first be 
            metaDataKeys: list of csvDict keys that should be put into
                          self.dbMetaData[row[1]]
            
        """
        for row in self.csvReader:
            if row[1] not in self.db: #row[1] is the chemical name
                self.db.update({row[1]:{}})
                self.db[row[1]].update({'name':row[1]})
                self.db[row[1]].update({'CAS':row[2]})
                self.db[row[1]].update({'formula':row[0]})
                self.dbMetaData.update({row[1]:{}})
                
            for key, val in csvDict.items():
                if len(val['colIdx']) == 1: #one col in row to process
                    try: #see the col can be converted to a float
                        self.db[row[1]].update({key:float(row[val['colIdx'][0]])})
                    except: #otherwise leave it alone (string)
                        self.db[row[1]].update({key:row[val['colIdx'][0]]})
                else: #multiple cols of floats in row to process
                    self.db[row[1]].update({key:[(float(row[idx])) for idx in val['colIdx']]})
                                
                if key not in self.dbMetaData[row[1]]: 
                    self.dbMetaData[row[1]].update({key:{}})
                for valKey,valVal in val.items():
                    if valKey in metaDataKeys:
                        self.dbMetaData[row[1]][key].update({valKey:valVal})
               
    def energyFormationGAS(self):
        """
        see doc strings for criticalProperty
        """
        csvFQPN = os.path.join(self.csvPath,'energyFormationGAS.csv')
        self.csvFile(csvFQPN)
        header = next(self.csvReader)
        header = [re.sub(self.refndSpace,'',h) for h in header]

        csvDict = {'eFGASState':
                       {'colIdx':[4],
                        'units':None,
                        'description':'energy Formation GAS State'},
                   'eFG298K':
                       {'colIdx':[5],
                        'units':self.reUnitSep.match(header[5]).groups()[1],
                        'description':'energy Formation, Gibbs at 298 K'},
                   'eFA298K':
                       {'colIdx':[6],
                        'units':self.reUnitSep.match(header[6]).groups()[1],
                        'description':'energy Formation, Helmholtz at 298 K'},
                   'eFGAS298K':
                       {'colIdx':[7],
                        'units':self.reUnitSep.match(header[7]).groups()[1],
                        'description':'energy Formation, Entropy at 298 K'}}
        #update csvDict for entries common to all top level keys
        for key in csvDict.keys():
            csvDict[key].update({'importDateTime':self.importDateTime,
                                 'dataSource':csvFQPN})

        metaDataKeys = ['description','dataSource','importDateTime','units']
        self.csv2db(csvDict,metaDataKeys)
        self.csvFileHandle.close() #close the file openend by csvFile

    def energyFormationHUS(self):
        """
        see doc strings for criticalProperty
        """
        csvFQPN = os.path.join(self.csvPath,'energyFormationHUS.csv')
        self.csvFile(csvFQPN)
        header = next(self.csvReader)
        header = [re.sub(self.refndSpace,'',h) for h in header]

        csvDict = {'eFHUSState':
                       {'colIdx':[4],
                        'units':None,
                        'description':'energy Formation HUS State'},
                   'eFH298K':
                       {'colIdx':[5],
                        'units':self.reUnitSep.match(header[5]).groups()[1],
                        'description':'energy Formation, Enthalpy at 298 K'},
                   'eFU298K':
                       {'colIdx':[6],
                        'units':self.reUnitSep.match(header[6]).groups()[1],
                        'description':'energy Formation, Internal at 298 K'},
                   'eFGAS298K':
                       {'colIdx':[7],
                        'units':self.reUnitSep.match(header[7]).groups()[1],
                        'description':'energy Formation, Entropy at 298 K'}}
        #update csvDict for entries common to all top level keys
        for key in csvDict.keys():
            csvDict[key].update({'importDateTime':self.importDateTime,
                                 'dataSource':csvFQPN})

        metaDataKeys = ['description','dataSource','importDateTime','units']
        self.csv2db(csvDict,metaDataKeys)
        self.csvFileHandle.close() #close the file openend by csvFile

    def enthalpyFusion(self):
        """
        see doc strings for criticalProperty
        """
        csvFQPN = os.path.join(self.csvPath,'enthalpyFusion.csv')
        self.csvFile(csvFQPN)
        header = next(self.csvReader)
        header = [re.sub(self.refndSpace,'',h) for h in header]

        csvDict = {'eFNFP':
                       {'colIdx':[5],
                        'units':self.reUnitSep.match(header[5]).groups()[1],
                        'description':'enthalpy Fusion at Normal Freezing Point'}}
        #update csvDict for entries common to all top level keys
        for key in csvDict.keys():
            csvDict[key].update({'importDateTime':self.importDateTime,
                                 'dataSource':csvFQPN})

        metaDataKeys = ['description','dataSource','importDateTime','units']
        self.csv2db(csvDict,metaDataKeys)
        self.csvFileHandle.close() #close the file openend by csvFile

    def enthalpyVaporization(self):
        """
        see doc strings for criticalProperty
        """
        csvFQPN = os.path.join(self.csvPath,'enthalpyVaporization.csv')
        self.csvFile(csvFQPN)
        header = next(self.csvReader)
        header = [re.sub(self.refndSpace,'',h) for h in header]

        csvDict = {'eVP':
                       {'colIdx':[3,4,5],
                        'units':'kJ/mol',
                        'description':'enthalpy Vaporization Parameters',
                        'eqnDoc':enthalpyVaporization.__doc__},
                   'eVTmin':
                       {'colIdx':[6],
                        'units':self.reUnitSep.match(header[6]).groups()[1],
                        'description':'enthalpy Vaporization Temperature Min'},
                    'eVTmax':
                       {'colIdx':[7],
                        'units':self.reUnitSep.match(header[7]).groups()[1],
                        'description':'enthalpy Vaporization Temperature Max'}}
        #update csvDict for entries common to all top level keys
        for key in csvDict.keys():
            csvDict[key].update({'importDateTime':self.importDateTime,
                                 'dataSource':csvFQPN})

        metaDataKeys = ['description','dataSource','importDateTime',
                        'units','eqnDoc']
        self.csv2db(csvDict,metaDataKeys)
        self.csvFileHandle.close() #close the file openend by csvFile
           
    def liquidDensity(self):
        """
        see doc strings for criticalProperty
        """
        csvFQPN = os.path.join(self.csvPath,'liquidDensity.csv')
        self.csvFile(csvFQPN)
        header = next(self.csvReader)
        header = [re.sub(self.refndSpace,'',h) for h in header]

        csvDict = {'lDP':
                       {'colIdx':[3,4,5,6],
                        'units':'g/cm^3',
                        'description':'liquid Denisty Parameters',
                        'eqnDoc':liquidDensity.__doc__},
                   'lDTmin':
                       {'colIdx':[7],
                        'units':self.reUnitSep.match(header[7]).groups()[1],
                        'description':'liquid Density Temperature Min'},
                    'lDTmax':
                       {'colIdx':[8],
                        'units':self.reUnitSep.match(header[8]).groups()[1],
                        'description':'liquid Density Temperature Max'}}
        #update csvDict for entries common to all top level keys
        for key in csvDict.keys():
            csvDict[key].update({'importDateTime':self.importDateTime,
                                 'dataSource':csvFQPN})

        metaDataKeys = ['description','dataSource','importDateTime',
                        'units','eqnDoc']
        self.csv2db(csvDict,metaDataKeys)
        self.csvFileHandle.close() #close the file openend by csvFile
                
    def liquidHeatCapacity(self):
        """
        see doc strings for criticalProperty
        """
        csvFQPN = os.path.join(self.csvPath,'liquidHeatCapacity.csv')
        self.csvFile(csvFQPN)
        header = next(self.csvReader)
        header = [re.sub(self.refndSpace,'',h) for h in header]

        csvDict = {'lHCP':
                       {'colIdx':[3,4,5,6],
                        'units':'J/mol/K',
                        'description':'liquid Heat Capacity Parameters',
                        'eqnDoc':liquidHeatCapacity.__doc__},
                   'lHCTmin':
                       {'colIdx':[7],
                        'units':self.reUnitSep.match(header[7]).groups()[1],
                        'description':'liquid Heat Capacity Temperature Min'},
                    'lHCTmax':
                       {'colIdx':[8],
                        'units':self.reUnitSep.match(header[8]).groups()[1],
                        'description':'liquid Heat Capacity Temperature Max'}}
        #update csvDict for entries common to all top level keys
        for key in csvDict.keys():
            csvDict[key].update({'importDateTime':self.importDateTime,
                                 'dataSource':csvFQPN})

        metaDataKeys = ['description','dataSource','importDateTime',
                        'units','eqnDoc']
        self.csv2db(csvDict,metaDataKeys)
        self.csvFileHandle.close() #close the file openend by csvFile

    def liquidThermalConductivity(self):
        """
        see doc strings for criticalProperty
        """
        csvFQPN = os.path.join(self.csvPath,'liquidThermalConductivity.csv')
        self.csvFile(csvFQPN)
        header = next(self.csvReader)
        header = [re.sub(self.refndSpace,'',h) for h in header]

        csvDict = {'lTCP':
                       {'colIdx':[3,4,5],
                        'units':'W/m/K',
                        'description':'liquid Thermal Conductivity Parameters',
                        'eqnDoc':liquidThermalConductivity.__doc__},
                   'lTCTmin':
                       {'colIdx':[6],
                        'units':self.reUnitSep.match(header[6]).groups()[1],
                        'description':'liquid Thermal Conductivity Temperature Min'},
                   'lTCTmax':
                       {'colIdx':[7],
                        'units':self.reUnitSep.match(header[7]).groups()[1],
                        'description':'liquid Thermal Conductivity Temperature Max'}}
        #update csvDict for entries common to all top level keys
        for key in csvDict.keys():
            csvDict[key].update({'importDateTime':self.importDateTime,
                                 'dataSource':csvFQPN})

        metaDataKeys = ['description','dataSource','importDateTime',
                        'units','eqnDoc']
        self.csv2db(csvDict,metaDataKeys)
        self.csvFileHandle.close() #close the file openend by csvFile
        
    def liquidViscosity(self):
        """
        see doc strings for criticalProperty
        """
        csvFQPN = os.path.join(self.csvPath,'liquidViscosity.csv')
        self.csvFile(csvFQPN)
        header = next(self.csvReader)
        header = [re.sub(self.refndSpace,'',h) for h in header]

        csvDict = {'lVP':
                       {'colIdx':[3,4,5,6],
                        'units':'cP',
                        'description':'liquid Viscosity Parameters',
                        'eqnDoc':liquidViscosity.__doc__},
                   'lVTmin':
                       {'colIdx':[7],
                        'units':self.reUnitSep.match(header[7]).groups()[1],
                        'description':'liquid Viscosity Temperature Min'},
                    'lVTmax':
                       {'colIdx':[8],
                        'units':self.reUnitSep.match(header[8]).groups()[1],
                        'description':'liquid Viscoisty Temperature Max'}}
        #update csvDict for entries common to all top level keys
        for key in csvDict.keys():
            csvDict[key].update({'importDateTime':self.importDateTime,
                                 'dataSource':csvFQPN})

        metaDataKeys = ['description','dataSource','importDateTime',
                        'units','eqnDoc']
        self.csv2db(csvDict,metaDataKeys)
        self.csvFileHandle.close() #close the file openend by csvFile

    def solidHeatCapacity(self):
        """
        see doc strings for criticalProperty
        """
        csvFQPN = os.path.join(self.csvPath,'solidHeatCapacity.csv')
        self.csvFile(csvFQPN)
        header = next(self.csvReader)
        header = [re.sub(self.refndSpace,'',h) for h in header]

        csvDict = {'sHCP':
                       {'colIdx':[3,4,5],
                        'units':'J/mol/K',
                        'description':'solid Heat Capacity Parameters',
                        'eqnDoc':solidHeatCapacity.__doc__},
                   'sHCTmin':
                       {'colIdx':[6],
                        'units':self.reUnitSep.match(header[6]).groups()[1],
                        'description':'solid Heat Capacity Temperature Min'},
                    'sHCTmax':
                       {'colIdx':[7],
                        'units':self.reUnitSep.match(header[7]).groups()[1],
                        'description':'solid Heat Capacity Temperature Max'}}
        #update csvDict for entries common to all top level keys
        for key in csvDict.keys():
            csvDict[key].update({'importDateTime':self.importDateTime,
                                 'dataSource':csvFQPN})

        metaDataKeys = ['description','dataSource','importDateTime',
                        'units','eqnDoc']
        self.csv2db(csvDict,metaDataKeys)
        self.csvFileHandle.close() #close the file openend by csvFile

    def surfaceTension(self):
        """
        see doc strings for criticalProperty
        """
        csvFQPN = os.path.join(self.csvPath,'surfaceTension.csv')
        self.csvFile(csvFQPN)
        header = next(self.csvReader)
        header = [re.sub(self.refndSpace,'',h) for h in header]

        csvDict = {'sTP':
                       {'colIdx':[3,4,5],
                        'units':'dyne/cm',
                        'description':'surface Tension Parameters',
                        'eqnDoc':surfaceTension.__doc__},
                   'sTTmin':
                       {'colIdx':[6],
                        'units':self.reUnitSep.match(header[6]).groups()[1],
                        'description':'surface Tension Temperature Min'},
                   'sTTmax':
                       {'colIdx':[7],
                        'units':self.reUnitSep.match(header[7]).groups()[1],
                        'description':'surface Tension Temperature Max'}}
        #update csvDict for entries common to all top level keys
        for key in csvDict.keys():
            csvDict[key].update({'importDateTime':self.importDateTime,
                                 'dataSource':csvFQPN})

        metaDataKeys = ['description','dataSource','importDateTime',
                        'units','eqnDoc']
        self.csv2db(csvDict,metaDataKeys)
        self.csvFileHandle.close() #close the file openend by csvFile
        
    def thermalConductivity(self):
        """
        see doc strings for criticalProperty
        """
        csvFQPN = os.path.join(self.csvPath,'thermalConductivity.csv')
        self.csvFile(csvFQPN)
        header = next(self.csvReader)
        header = [re.sub(self.refndSpace,'',h) for h in header]

        csvDict = {'tCP':
                       {'colIdx':[3,4,5],
                        'units':'W/m/K',
                        'description':'thermal Conductivity Parameters',
                        'eqnDoc':thermalConductivity.__doc__},
                   'tCTmin':
                       {'colIdx':[6],
                        'units':self.reUnitSep.match(header[6]).groups()[1],
                        'description':'thermal Conductivity Temperature Min'},
                   'tCTmax':
                       {'colIdx':[7],
                        'units':self.reUnitSep.match(header[7]).groups()[1],
                        'description':'thermal Conductivity Temperature Max'},
                   'tCState':
                       {'colIdx':[8],
                        'units':None,
                        'description':'thermal Conductivity State'}}
        #update csvDict for entries common to all top level keys
        for key in csvDict.keys():
            csvDict[key].update({'importDateTime':self.importDateTime,
                                 'dataSource':csvFQPN})

        metaDataKeys = ['description','dataSource','importDateTime',
                        'units','eqnDoc']
        self.csv2db(csvDict,metaDataKeys)
        self.csvFileHandle.close() #close the file openend by csvFile

    def vaporEnthalpyFormation(self):
        """
        see doc strings for criticalProperty
        """
        csvFQPN = os.path.join(self.csvPath,'vaporEnthalpyFormation.csv')
        self.csvFile(csvFQPN)
        header = next(self.csvReader)
        header = [re.sub(self.refndSpace,'',h) for h in header]

        csvDict = {'vEFP':
                       {'colIdx':[3,4,5],
                        'units':'kJ/mol',
                        'description':'vapor Enthalpy Formation Parameters',
                        'eqnDoc':vaporEnthalpyFormation.__doc__},
                   'vEFTmin':
                       {'colIdx':[6],
                        'units':self.reUnitSep.match(header[6]).groups()[1],
                        'description':'vapor Enthalpy Formation Temperature Min'},
                   'vEFTmax':
                       {'colIdx':[7],
                        'units':self.reUnitSep.match(header[7]).groups()[1],
                        'description':'vapor Enthalpy Formation Temperature Max'}}
        #update csvDict for entries common to all top level keys
        for key in csvDict.keys():
            csvDict[key].update({'importDateTime':self.importDateTime,
                                 'dataSource':csvFQPN})

        metaDataKeys = ['description','dataSource','importDateTime',
                        'units','eqnDoc']
        self.csv2db(csvDict,metaDataKeys)
        self.csvFileHandle.close() #close the file openend by csvFile

    def vaporGibbsEnergyFormation(self):
        """
        see doc strings for criticalProperty
        """
        csvFQPN = os.path.join(self.csvPath,'vaporGibbsEnergyFormation.csv')
        self.csvFile(csvFQPN)
        header = next(self.csvReader)
        header = [re.sub(self.refndSpace,'',h) for h in header]

        csvDict = {'vGEFP':
                       {'colIdx':[3,4,5],
                        'units':'kJ/mol',
                        'description':'vapor Gibbs Energy Formation Parameters',
                        'eqnDoc':vaporGibbsEnergyFormation.__doc__},
                   'vGEFTmin':
                       {'colIdx':[6],
                        'units':self.reUnitSep.match(header[6]).groups()[1],
                        'description':'vapor Gibbs Energy Formation Temperature Min'},
                   'vGEFTmax':
                       {'colIdx':[7],
                        'units':self.reUnitSep.match(header[7]).groups()[1],
                        'description':'vapor Gibbs Energy Formation Temperature Max'}}
        #update csvDict for entries common to all top level keys
        for key in csvDict.keys():
            csvDict[key].update({'importDateTime':self.importDateTime,
                                 'dataSource':csvFQPN})

        metaDataKeys = ['description','dataSource','importDateTime',
                        'units','eqnDoc']
        self.csv2db(csvDict,metaDataKeys)
        self.csvFileHandle.close() #close the file openend by csvFile
        
    def vaporHeatCapacity(self):
        """
        see doc strings for criticalProperty
        """
        csvFQPN = os.path.join(self.csvPath,'vaporHeatCapacity.csv')
        self.csvFile(csvFQPN)
        header = next(self.csvReader)
        header = [re.sub(self.refndSpace,'',h) for h in header]

        csvDict = {'vHCP':
                       {'colIdx':[3,4,5,6,7],
                        'units':'J/mol/K',
                        'description':'vapor Heat Capacity Parameters',
                        'eqnDoc':vaporHeatCapacity.__doc__},
                   'vHCTmin':
                       {'colIdx':[8],
                        'units':self.reUnitSep.match(header[8]).groups()[1],
                        'description':'vapor Heat Capacity Temperature Min'},
                    'vHCTmax':
                       {'colIdx':[9],
                        'units':self.reUnitSep.match(header[9]).groups()[1],
                        'description':'vapor Heat Capacity Temperature Max'}}
        #update csvDict for entries common to all top level keys
        for key in csvDict.keys():
            csvDict[key].update({'importDateTime':self.importDateTime,
                                 'dataSource':csvFQPN})

        metaDataKeys = ['description','dataSource','importDateTime',
                        'units','eqnDoc']
        self.csv2db(csvDict,metaDataKeys)
        self.csvFileHandle.close() #close the file openend by csvFile

    def vaporPressure(self):
        """
        see doc strings for criticalProperty
        """
        csvFQPN = os.path.join(self.csvPath,'vaporPressure.csv')
        self.csvFile(csvFQPN)
        header = next(self.csvReader)
        header = [re.sub(self.refndSpace,'',h) for h in header]

        csvDict = {'vPP':
                       {'colIdx':[3,4,5,6,7],
                        'units':'mmHg',
                        'description':'vapor Pressure Parameters',
                        'eqnDoc':vaporPressure.__doc__},
                   'vPTmin':
                       {'colIdx':[8],
                        'units':self.reUnitSep.match(header[8]).groups()[1],
                        'description':'vapor Pressure Temperature Min'},
                    'vPTmax':
                       {'colIdx':[9],
                        'units':self.reUnitSep.match(header[9]).groups()[1],
                        'description':'vapor Pressure Temperature Max'}}
        #update csvDict for entries common to all top level keys
        for key in csvDict.keys():
            csvDict[key].update({'importDateTime':self.importDateTime,
                                 'dataSource':csvFQPN})

        metaDataKeys = ['description','dataSource','importDateTime',
                        'units','eqnDoc']
        self.csv2db(csvDict,metaDataKeys)
        self.csvFileHandle.close() #close the file openend by csvFile
                
    def vaporThermalConductivity(self):
        """
        see doc strings for criticalProperty
        """
        csvFQPN = os.path.join(self.csvPath,'vaporThermalConductivity.csv')
        self.csvFile(csvFQPN)
        header = next(self.csvReader)
        header = [re.sub(self.refndSpace,'',h) for h in header]

        csvDict = {'vTCP':
                       {'colIdx':[3,4,5],
                        'units':'W/m/K',
                        'description':'organic Vapor Thermal Conductivity Parameters',
                        'eqnDoc':vaporThermalConductivity.__doc__},
                   'vTCTmin':
                       {'colIdx':[6],
                        'units':self.reUnitSep.match(header[6]).groups()[1],
                        'description':'organic Vapor Thermal Conductivity Temperature Min'},
                   'vTCTmax':
                       {'colIdx':[7],
                        'units':self.reUnitSep.match(header[7]).groups()[1],
                        'description':'organic Vapor Thermal Conductivity Temperature Max'}}
        #update csvDict for entries common to all top level keys
        for key in csvDict.keys():
            csvDict[key].update({'importDateTime':self.importDateTime,
                                 'dataSource':csvFQPN})

        metaDataKeys = ['description','dataSource','importDateTime',
                        'units','eqnDoc']
        self.csv2db(csvDict,metaDataKeys)
        self.csvFileHandle.close() #close the file openend by csvFile

    def vaporViscosity(self):
        """
        see doc strings for criticalProperty
        """
        csvFQPN = os.path.join(self.csvPath,'vaporViscosity.csv')
        self.csvFile(csvFQPN)
        header = next(self.csvReader)
        header = [re.sub(self.refndSpace,'',h) for h in header]

        csvDict = {'vVP':
                       {'colIdx':[3,4,5],
                        'units':'micropoise',
                        'description':'vapor Viscosity Parameters',
                        'eqnDoc':vaporViscosity.__doc__},
                   'vVTmin':
                       {'colIdx':[6],
                        'units':self.reUnitSep.match(header[6]).groups()[1],
                        'description':'vapor Viscosity Temperature Min'},
                   'vVTmax':
                       {'colIdx':[7],
                        'units':self.reUnitSep.match(header[7]).groups()[1],
                        'description':'vapor Viscosity Temperature Max'}}
        #update csvDict for entries common to all top level keys
        for key in csvDict.keys():
            csvDict[key].update({'importDateTime':self.importDateTime,
                                 'dataSource':csvFQPN})

        metaDataKeys = ['description','dataSource','importDateTime',
                        'units','eqnDoc']
        self.csv2db(csvDict,metaDataKeys)
        self.csvFileHandle.close() #close the file openend by csvFile
        
    def waterSolubility(self):
        """
        see doc strings for criticalProperty
        """
        csvFQPN = os.path.join(self.csvPath,'waterSolubility.csv')
        self.csvFile(csvFQPN)
        header = next(self.csvReader)
        header = [re.sub(self.refndSpace,'',h) for h in header]

        csvDict = {'wSP':
                       {'colIdx':[5,6,7],
                        'units':'ppm',
                        'description':'water Solubility Parameters',
                        'eqnDoc':waterSolubility.__doc__}}
        #update csvDict for entries common to all top level keys
        for key in csvDict.keys():
            csvDict[key].update({'importDateTime':self.importDateTime,
                                 'dataSource':csvFQPN})

        metaDataKeys = ['description','dataSource','importDateTime',
                        'units','eqnDoc']
        self.csv2db(csvDict,metaDataKeys)
        self.csvFileHandle.close() #close the file openend by csvFile
           
    def saltWaterSolubility(self):
        """
        see doc strings for criticalProperty
        """
        csvFQPN = os.path.join(self.csvPath,'saltWaterSolubility.csv')
        self.csvFile(csvFQPN)
        header = next(self.csvReader)
        header = [re.sub(self.refndSpace,'',h) for h in header]

        csvDict = {'sWSP':
                       {'colIdx':[5,6,7],
                        'units':'ppm',
                        'description':'salt Water Solubility Parameters',
                        'eqnDoc':saltWaterSolubility.__doc__}}
        #update csvDict for entries common to all top level keys
        for key in csvDict.keys():
            csvDict[key].update({'importDateTime':self.importDateTime,
                                 'dataSource':csvFQPN})

        metaDataKeys = ['description','dataSource','importDateTime',
                        'units','eqnDoc']
        self.csv2db(csvDict,metaDataKeys)
        self.csvFileHandle.close() #close the file openend by csvFile
                   
if __name__ == "__main__":
    ycphi = YCPH()
    [ydb, ydbmd] = ycphi.importYaws()
    
    saveFileName='yaws4mpdb-' + strftime("%Y%m%d-%H%M") + '.pickle'
    saveFileName=os.path.join('yaws',saveFileName)
    dbFileName = dbSave(ydb,ydbmd,fileName=saveFileName)
    print("Created yaws mpdb pickle: %s" % dbFileName)
