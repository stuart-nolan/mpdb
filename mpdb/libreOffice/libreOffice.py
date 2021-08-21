"""
mpdb: Material Property Data Base as python module
libreOffice.py: utility module for interacting with libreoffice

Revision Date: 2021.08.18

SPDX-License-Identifier: BSD-2-Clause
Copyright (c) 2021 Stuart Nolan. All rights reserved.
"""
import uno
import os
import sys
import platform
from subprocess import Popen, PIPE
from time import sleep
from mpdb.retry import retry

class libreOfficeInterface():
    """
    Convenience functions for interfacing with libreoffice
    
    start libreoffice like: 
        soffice --calc --accept="socket,host=localhost,port=2002;urp;StarOffice.ServiceManager"
    
    Ref:
        http://christopher5106.github.io/office/2015/12/06/openoffice-libreoffice-automate-your-office-tasks-with-python-macros.html
    """
    def __init__(self,documentName=None):
        # get the uno component context from the PyUNO runtime
        self.localContext = uno.getComponentContext()
        # create the UnoUrlResolver
        self.resolver = self.localContext.ServiceManager.createInstanceWithContext(
        				"com.sun.star.bridge.UnoUrlResolver", self.localContext )
        # connect to the running office
        self.document = self.initActiveDoc()
        #self.document.Title
        self.activeSheet = self.getSheet()
        #self.activeSheet.Name
    
    @retry(tries=3)
    def initActiveDoc(self):
        """
        return the currently active libre office document
        
        use this function with @retry just after starting soffice to give 
        soffice time to load
        """
        self.ctx = self.resolver.resolve( "uno:socket,host=localhost,port=2002;urp;StarOffice.ComponentContext")
        self.smgr = self.ctx.ServiceManager
        # get the central desktop object
        self.desktop = self.smgr.createInstanceWithContext( "com.sun.star.frame.Desktop",self.ctx)
        # access the current document but need to wait a few seconds for 
        #soffice to load...

        activeDoc = self.desktop.getCurrentComponent()
        if activeDoc is None: 
            raise Exception("No Active Document")
        else:
            return activeDoc

    def getDocument(self,docTitle='__activeDoc__'):
        """
        return a libre office document object

        Parameters:
            docTitle, string name of an open libre office document
        
        Notes:
            if docTitle is not in ?? then return active document
            if no docTitle or active document, msg user
            [doc.Title for doc in self.desktop.getComponents().createEnumeration()]
        """
        docTitles = [doc.Title for doc in self.desktop.getComponents()]
        if activeDoc is None: 
            raise Exception("No Active Document")
        else:
            return activeDoc

    def getSheet(self,sheetName='__activeSheet__'):
        """
        return a libre office Sheet object

        Parameters:
            sheetName, string name of a sheet in an open libre office document
        
        Notes:
            if sheetName is not in document.Sheets then return active sheet            
        """
        try:
            if self.document.Sheets.hasByName(sheetName):
                return self.document.Sheets.getByName(sheetName)
            else:
                return self.document.CurrentController.ActiveSheet
        except AttributeError:
            return None
    
    def setName(self,cellRangeName,name,sheet):
        """
        set a named range
        
        Parameters:
            cellRangeName, string value of absolute cell range, i.e. '$C$2'
            name, string to name new cell range
            sheet, uno libre office Sheet object (see getSheet)
        
        Notes:
            if name exists then nothing is done
        """
        try:
            if not self.document.NamedRanges.hasByName(name):
                oCellAddress = sheet.getCellRangeByName(cellRangeName).getCellAddress()
                self.document.NamedRanges.addNewByName(name,cellRangeName,oCellAddress,0)
        except AttributeError:
            return None
    
    def listNames(self):
        """
        list named ranges
        
        Returns:
            tuple of strings of named ranges
        Notes:
            
        """
        try:
            return self.document.NamedRanges.getElementNames()
        except AttributeError:
            return ()

    def getNamedRange(self,name):
        """
        Returns:
            libre office named range object

        Notes:
            self.setName("$C$3","fw_water",self.getSheet())
            namedRange = getNamedRange("fw_water")
            namedRange.getContent() returns '$C$3'
            namedRange.ReferredCells.Value = 18 sets cell C3 to 18
            
        """
        try:
            if self.document.NamedRanges.hasByName(name):
                return self.document.NamedRanges.getByName(name)
        except AttributeError:
            return None
        
def libreOfficeCalc():
    """
    Ref: https://stackoverflow.com/a/13256908/4279
    """
    # set system/version dependent "start_new_session" analogs
    kwargs = {}
    if platform.system() == 'Windows':
        # from msdn [1]
        CREATE_NEW_PROCESS_GROUP = 0x00000200  # note: could get it from subprocess
        DETACHED_PROCESS = 0x00000008          # 0x8 | 0x200 == 0x208
        kwargs.update(creationflags=DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP)  
    elif sys.version_info < (3, 2):  # assume posix
        kwargs.update(preexec_fn=os.setsid)
    else:  # Python 3.2+ and Unix
        kwargs.update(close_fds=True)

    p = Popen(["soffice", "--calc",
               "--accept=socket,host=localhost,port=2002;urp;StarOffice.ServiceManager"], stdin=None, stdout=None, stderr=None, **kwargs)
    assert not p.poll()
    #return os.system('soffice --calc --accept="socket,host=localhost,port=2002;urp;StarOffice.ServiceManager"&')

if __name__ == "__main__":
    libreOfficeCalc()
    lOI = libreOfficeInterface()
    """
    python scripts run from inside lo refer to XSCRIPTCONTEXT, that object can
    be created remotely (i.e. in this script) via:

    from pythonscript import ScriptContext
    XSCRIPTCONTEXT = ScriptContext(lOI.ctx, None, None)
    doc = XSCRIPTCONTEXT.getDocument()

    """
    
    #activeSheet = lOI.getSheet()
    
    # access cell C4
    cell1 = lOI.activeSheet.getCellRangeByName("C4")
    
    # set text inside
    cell1.String = "Hello world"
    
    # other example with a value
    cell2 = lOI.activeSheet.getCellRangeByName("E6")
    cell2.Value = cell2.Value + 1
    cell3 = lOI.activeSheet.getCellRangeByName("b2")
    cell3.Formula = "=E6"
    # access cell C4 on Sheet2
    cell4 = lOI.getSheet("Sheet1").getCellRangeByName("C5")
    
    # set text inside
    cell4.String = "Hello Sheet2"
