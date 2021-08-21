#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mpdb: Material Property Data Base as python module
mpdb-gui.py: material property data base gui

Author: Stuart Nolan
Revision Date: 2017/08/30
Copyright (c) 2006 Stuart Nolan

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

  Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
AS IS AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

TODO:
"""

import os
import sys
#import inspect
import json
from time import strftime, sleep
import argparse
import platform
from appdirs import *
from PyQt5 import QtCore, QtGui, QtWidgets

from mpdb.periodicTable import fW
from mpdb.units import *
from mpdb.pcmpEQ import *
from mpdb.dbUtil import dbFind, dbKeys, dbLoad, dbSave

location = os.path.dirname(os.path.realpath(os.path.abspath(__file__)))

from mpdbui import Ui_MainWindow
    
def setTrace():
    """
    Set a tracepoint in the Python debugger that works with Qt
    from PyQt4.QtCore import pyqtRemoveInputHook
    when done debugging, call QtCore.pyqtRestoreInputHook() from pdb
    may have to use c several times
    """

    from PyQt5.QtCore import pyqtRemoveInputHook

    from pdb import set_trace
    pyqtRemoveInputHook()
    set_trace()

def removeTrace():
    """
    """
    from PyQt5.QtCore import pyqtRestoreInputHook
    pyqtRestoreInputHook()
    
class mpdbUI(Ui_MainWindow):
    def __init__(self):
        """
        Workflow:
            edit & save form from Qt designer
            from bash shell:
                pyuic5 -x mpdb.ui -o mpdbui.py
                python mpdb-gui.py -p bob
        
        Parameters:
            opts, argparse dictionary.  Argparse is used in main() to create 
                  this dictionary from the command line.
        """
        Ui_MainWindow.__init__(self)
        self.__version__ = "0.1"
        self.confFileName = "mpdbrc.json"
        self.conf = {}
        self.confDir = user_data_dir("mpdb")
        
    def cmdLineArgs(self):
        """
        """
        self.ignorUserConf = self.opts['ignorUserConf']
        
    def loadConf(self):
        """
        """
        with open(self.confFQPN, 'r') as infile:
            self.conf = json.load(infile)
    
    def defaultConf(self):
        """
        """
        self.conf.update({'version':self.__version__})
    
    def saveConf(self):
        """
        """
        with open(self.confFQPN, 'w') as outfile:
            json.dump(self.conf, outfile)
        
    def start(self,opts):
        """
        Call after __init__ and setupUI to process command line arguments, 
        read the user configuration file (or create it if it does not exist), 
        and connect signals to slot functions defined below
        """
        self.opts = opts
        self.cmdLineArgs()
        self.confFQPN = os.path.join(self.confDir,self.confFileName)
        if not self.ignorUserConf:
            if not os.path.exists(self.confDir):
                os.makedirs(self.confDir)
            if not os.path.isfile(self.confFQPN):
                self.defaultConf()
                self.saveConf()
            self.loadConf()
        else:
            self.defaultConf()

        self.dbModel = QtGui.QStandardItemModel()
        self.dbModel.setHorizontalHeaderLabels(['Name','Forumla','CAS'])
        self.proxyModel = QtCore.QSortFilterProxyModel(self.treeView)
        self.proxyModel.setSourceModel(self.dbModel)
        self.proxyModel.setFilterKeyColumn(-1) #-1 = all columns
        self.treeView.setSortingEnabled(True)
        self.treeView.setModel(self.proxyModel)
        self.lineEdit.textChanged.connect(self.onTextChanged)    
        self.treeView.clicked.connect(self.selectDBItem)

        self.dataModel = QtGui.QStandardItemModel()
        self.proxyModel_2 = QtCore.QSortFilterProxyModel(self.treeView_2)
        self.proxyModel_2.setSourceModel(self.dataModel)
        #self.proxyModel_2.setFilterKeyColumn(-1) #-1 = all columns
        self.treeView_2.setSortingEnabled(True)
        self.treeView_2.setModel(self.proxyModel_2)
        self.lineEdit_2.textChanged.connect(self.onTextChanged_2)    
        self.dataModel.setHorizontalHeaderLabels(['ID','Value','Units'])
        self.action_Quit.triggered.connect(self.slot_Quit)
        self.action_Open.triggered.connect(self.slot_Open)
        
    def slot_Quit(self):
        """
        """
        try:
            self.saveConf()
        except:
            print("Error on mpdb GUI exit.  Possible data loss.")
        sys.exit(0)
        
    def slot_Open(self):
        """
        """
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        self.fileName, _ = QtWidgets.QFileDialog.getOpenFileName(None,"Open Database Pickle", "","All Files (*);;Pickle Files (*.pickle)", options=options)
        #self.statusbar.showMessage(self.fileName,2000)
        self.openDB()
            
    def openDB(self):
        """
        """
        (self.db, self.dbMetaData) = dbLoad(fileName=self.fileName,interactive=False)
        self.dbNames, self.dbCAS, self.dbFormula = zip(*dbKeys(self.db))
        self.groupBox.setTitle(self.fileName)

        for idx,name in enumerate(self.dbNames):
            name_item = QtGui.QStandardItem(name)
            self.dbModel.appendRow([name_item,
                                    QtGui.QStandardItem(self.dbFormula[idx]),
                                    QtGui.QStandardItem(self.dbCAS[idx])])
    
    def onTextChanged(self, text):
        self.proxyModel.setFilterRegExp(text)

    def onTextChanged_2(self, text):
        self.proxyModel_2.setFilterRegExp(text)

    def selectDBItem(self, index):
        ix = self.proxyModel.mapToSource(index)
        parent = self.dbModel.itemFromIndex(ix)
        name = parent.text()
        self.dataModel.clear()
        self.dataModel.setHorizontalHeaderLabels(['ID','Value','Units'])

        if name not in self.db:
            self.statusbar.showMessage("%s not found, try clicking on a name"
                                       % name,5000)
            return

        for key,val in self.db[name].items():
            key_item = QtGui.QStandardItem(key)
            val_item = QtGui.QStandardItem(str(val))
            units = ""
            if key in self.dbMetaData[name]:
                if 'units' in self.dbMetaData[name][key]:
                    units = self.dbMetaData[name][key]['units']

                for dbMDkey,dbMDval in self.dbMetaData[name][key].items():
                    dbMDkey_item = QtGui.QStandardItem(dbMDkey)
                    dbMDval_item = QtGui.QStandardItem(str(dbMDval))
                    key_item.appendRow([dbMDkey_item,dbMDval_item])

            units_item = QtGui.QStandardItem(units)
            self.dataModel.appendRow([key_item,val_item,units_item])

    
    def closeEvent(self):
        """
        """
        self.slot_Quit()
        sys.exit(0)

def main():
    """
    command line inputs for mpdb-gui  
    """
    ui = mpdbUI()
    p = argparse.ArgumentParser(description='GUI for mpdb')
    p.add_argument('-i', '--ignorUserConf', action='store_true',
                   help='print log info to stdout',default=False)    
    p.add_argument('--version', action='version', version='%(prog)s Version: ' +
                                                               ui.__version__)
    p.add_argument('-t', '--testVal', action='store', default='',
                   metavar='cmd line arg test value',
                   help='used for development and testing')
   
    args = p.parse_args()

    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui.setupUi(MainWindow)
    ui.start(vars(args))
    MainWindow.show()
    app.aboutToQuit.connect(ui.closeEvent)
    sys.exit(app.exec_())
    
if __name__ == '__main__':
    """
    TODO:
    """
    main()
