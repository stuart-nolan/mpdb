"""
mpdb: Material Property Data Base

Revision Date: 2021.08.21

SPDX-License-Identifier: BSD-2-Clause
Copyright (c) 2021 Stuart Nolan. All rights reserved.
"""
import json
import os
import appdirs
from pkgutil import extend_path

__path__ = extend_path(__path__, __name__)
location = os.path.dirname(os.path.realpath(os.path.abspath(__file__)))
__version__ = "0.1a0"
__all__ = ["costIndexes",
           "eq/misc_python",
           "eq/yaws_python",
           "gui/mpdbTk",
           "periodicTable",
           "units",
           "utils"]
__data__ =["data/slopImport",
           "data/yawsImport",
           "data/rSanderImport"]
__cython__=["eq/yaws_cython",
            "eq/misc_cython"]

prefFileName = "%s_rc.json" % __name__
prefDir = appdirs.user_data_dir(__name__)
prefFQPN = os.path.join(prefDir,prefFileName)

def defaultPref():
    return {'version':__version__,
            'lastOpenDir':'',
            'lastOpenFQPNs':[],
            'dataDir':'data',
            'costIndexFile':os.path.join('data',
                                         'costIndexes',
                                         'costIndexes.csv')}

pref = defaultPref()

def savePref(pref=pref,prefFQPN=prefFQPN):
    with open(prefFQPN, 'w') as outfile:
        json.dump(pref, outfile)

def loadPref(prefFQPN):
    if not os.path.exists(prefDir):
        os.makedirs(prefDir)
    if not os.path.isfile(prefFQPN):
        savePref()
    with open(prefFQPN, 'r') as infile:
        return json.load(infile)

pref = loadPref(prefFQPN)
