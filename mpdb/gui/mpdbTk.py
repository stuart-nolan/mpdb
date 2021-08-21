"""
mpdb: Material Property Data Base as python module
mpdbTk.py: material property data base tkinter gui

Revision Date: 2021.08.15

SPDX-License-Identifier: BSD-2-Clause
Copyright (c) 2021 Stuart Nolan. All rights reserved.

TODO:
context menu for data items
export data as csv
create named cells in libreOffice and/or excel workbooks?
export equations for libreOffice and/or excel (vba?)
data/eqn display
minimal help
"""

# BEGIN python standard library imports
import os
import sys
if sys.version_info[0] != 3 or sys.version_info[1] < 5:
    print("Python version >= 3.5 is required.")
    sys.exit(1)
sys.path.insert(0, "/home/ul/local/src/mpdb/")
import io
import csv
import pickle
from time import strftime, sleep
import tkinter as tk
from tkinter import ttk
from tkinter.filedialog import askopenfilename
#from pdb import set_trace
# set_trace()
# END python standard library imports

# BEGIN python foreign library imports
import argparse
# END python foreign library imports

# BEGIN app library imports
import mpdb
from mpdb.gui.tkMenu import tkMenubar
from mpdb.gui.tkStatusBar import StatusBar
from mpdb.gui.ttkNotebook import Notebook
from mpdb.gui.ttkFilterTree import FilterTree
from mpdb.periodicTable import fW
from mpdb.units import *
from mpdb.eq.yaws_python import *
# END app library imports

def isnotebook():
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return True   # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return False  # Terminal running IPython
        else:
            return False  # Other type (?)
    except NameError:
        return False      # Probably standard Python interpreter

class GUI(ttk.Frame):
    def __init__(self, parent, *args, **kwargs):
        """
        Init the app
        """
        ttk.Frame.__init__(self, parent, *args, **kwargs)
        # application parameters and prefs
        self.__version__ = mpdb.__version__
        self.appName = mpdb.__name__
        self.maxOpenFQPNs = 10
        self.db={}
        self.dbMetaData={}
        
        # main window title
        self.master.title(self.appName)
        
        # menubar & menus
        self.menubar = tkMenubar(parent,relief=tk.FLAT)
        self.master.config(menu=self.menubar)
        
        # menu commands (see tkMenubar.py for menu names and entry names)
        self.menubar.menus['File'].entryconfig('Open',command=self.open)
        self.menubar.menus['File'].entryconfig('Save',command=self.savePref)
        self.menubar.menus['File'].entryconfig('Quit',command=self.quit)
        self.menubar.menus['Edit'].entryconfig('Copy',command=self.copy)
        self.menubar.menus['Tools'].entryconfig('Clear Open Recent',
                                                command=self.clearRecents)
        
        # keyboard accelerators
        self.bind_all("<Control-o>", self.open)
        self.bind_all("<Control-q>", self.quit)
        self.bind_all("<Control-c>", self.copy)
        
        # status bar
        self.statusBar = StatusBar(parent)
        self.statusBar.set("%s","status bar initialized")
        self.statusBar.pack(side=tk.BOTTOM, fill=tk.X)
        
        # notebook
        self.tabDefs = [{"text":"Browse"},
                        {"text":"Prefs"}]
        self.nb = Notebook(parent,self.tabDefs)
        self.nb.pack(fill=tk.BOTH,expand=True)
        
        # filter tree view
        self.filterTree = FilterTree(self.nb.tabFrames[self.nb.tfids["Browse"]])
        self.filterTree.pack(fill='both', expand=True)
        
        # main window size & position
        self.master.withdraw()
        
        # screen width and height
        ws = parent.winfo_screenwidth() # width of the screen
        hs = parent.winfo_screenheight() # height of the screen
        wf=0.5
        hf=0.9
        # x pixels from right side, y pixels from top
        x = (ws)*(1 - wf)
        y = 0
        parent.geometry('%dx%d+%d+%d' % (ws*wf, hs*hf, x, y))
        parent.protocol("WM_DELETE_WINDOW", self.quit)

    def clearRecents(self):
        self.menubar.clearRecents()
        mpdb.pref['lastOpenFQPNs'] = []
        
    def cmdLineArgs(self):
        """
        """
        self.ignoreUserPref = self.opts['ignoreUserPref']
        self.showSplash = self.opts['showSplash']

    def copy(self, *args):
        self.clipboard_clear()
        valColMap = {'#1': 0, '#2': 1}
        if args:
            event = args[0]
            col = self.filterTree.tv.identify_column(event.x)
        else:
            col='#0'
            
        item = self.filterTree.tv.item(self.filterTree.tv.focus())
        if col == '#0':
            copyVal = item['text']
        else:
            copyVal = item['values'][valColMap[col]]
        self.clipboard_append(copyVal)
        
        """
        items = self.filterTree.tv.selection()
        cpItems = ''
        for item in items:
            if cpItems == '':
                cpItems = self.filterTree.tv.item(item)['text']
            else:
                cpItems = cpItems + '\r\n' + \
                          self.filterTree.tv.item(item)['text']

            if any(nl in self.filterTree.tv.item(item)['values'][0]
                   for nl in ['\n','\r\n']):
                sep = '\r\n'
            else:
                sep = ','
            
            cpItems = cpItems + \
                      sep + \
                      ','.join(map(str,self.filterTree.tv.item(item)['values']))
        """
        
    def csvStr(self,item):
        """
        returns:
         - a coma separated string from items in a list
        """
        if hasattr(item,'__iter__'):
            if any(isinstance(i,list) for i in item):
                out = io.StringIO()
                csvw = csv.writer(out)
                csvw.writerows(item)
                return out.getvalue()
            elif isinstance(item,list):
                out = io.StringIO()
                csvw = csv.writer(out,lineterminator="")
                csvw.writerow(item)
                return out.getvalue()
        return item

    def open(self, *args):
        """
        """
        if 'lastOpenDir' not in mpdb.pref.keys():
            lastDir=os.getcwd()
        elif os.path.isdir(mpdb.pref['lastOpenDir']):
            lastDir = mpdb.pref['lastOpenDir']
        else:
            lastDir=os.getcwd()
            
        FILEOPENOPTIONS = dict(defaultextension='.pickle',
                               initialdir=lastDir,
                               filetypes=[('pickle file','*.pickle'),
                                          ('All files','*.*')])
        fileName=askopenfilename(parent=self, **FILEOPENOPTIONS)
        
        mpdb.pref.update({'lastOpenDir':os.path.split(fileName)[0]})
        self.openDB(fileName)
            
    def openDB(self, fileName):
        """
        Opens a material database pickle file with the FQPN fileName,
        populates, the filter tree view, updates mpdb.pref['lastOpenFQPNs'], 
        and updates the 'Open Recent' menu.
        
        fileName pickle contains a two item list where each list item is 
        a dictionary.

        """
        try:
            with open(fileName, 'rb') as f:
                db = pickle.load(f)

            if type(db) in [list, tuple]:
                if type(db[0]) is dict and type(db[1]) is dict:
                    self.statusBar.set("Loading %s",fileName)
                    (self.db, self.dbMetaData) = db
                    tk.Tk.update(self)
                    self.populateFTV()
                    self.master.title(self.appName+": "+
                                      os.path.split(fileName)[1])
                    if 'lastOpenFQPNs' not in mpdb.pref.keys():
                        mpdb.pref.update({'lastOpenFQPNs':[fileName]})
                    elif len(mpdb.pref['lastOpenFQPNs']) > self.maxOpenFQPNs:
                        if fileName not in mpdb.pref['lastOpenFQPNs']:
                            mpdb.pref['lastOpenFQPNs'].pop(0)
                            mpdb.pref['lastOpenFQPNs'].append(fileName)
                    else:
                        if fileName not in mpdb.pref['lastOpenFQPNs']:
                            mpdb.pref['lastOpenFQPNs'].append(fileName)
                            self.menubar.updateRecents(mpdb.pref['lastOpenFQPNs'],
                                               getattr(self,'openDB'))
                else:
                    self.statusBar.set("Invalid file %s: ",fileName)
            else:
                self.statusBar.set("Invalid file %s: ",fileName)

        except IOError:
            self.statusBar.set("Error trying to open %s",fileName)

        except Exception as e:
            print(e)
            self.statusBar.set("Exception opening %s",fileName)
        
    def populateFTV(self):
        self.filterTree.tv.heading('#0',text='Name', command=lambda _col='#0':
                                   self.filterTree.tvSortColumn(_col, False))
        self.columnDef = {'for':'Formula or Value','cas':'CAS or Units'}
        self.filterTree.tv.configure(columns=[k for k in self.columnDef.keys()])
        for key,val in self.columnDef.items():
            #self.filterTree.tv.column(key,stretch=0)
            self.filterTree.tv.heading(key,
                                       text=key,
                                       command=lambda
                                       colID = key:
                                       self.filterTree.tvSortColumn(colID,
                                                                    False))
            self.filterTree.tv.heading(key, text=val)
        self.filterTree.tv.column('#0',
                                  width=int(self.filterTree.tv.winfo_width()/3))
        #delete preveious entries
        self.filterTree.tv.delete(*self.filterTree.tv.get_children())
        #insert new data from db
        for key,val in self.db.items():
            if 'CAS' not in val.keys():
                vals=(val['formula'],"")
            else:
                vals=(val['formula'],val['CAS'])
            vID = self.filterTree.tv.insert('', 'end', text=key, values=vals)
            for vKey,vVal in val.items():
                vVal=self.csvStr(vVal)
                vals=(vVal,"")
                if vKey in self.dbMetaData[key].keys():
                    if 'units' in self.dbMetaData[key][vKey].keys():
                        vals=(vVal,self.dbMetaData[key][vKey]['units'])
                vValID = self.filterTree.tv.insert(vID,
                                                   'end',
                                                   text=vKey,
                                                   values=vals)
                if vKey in self.dbMetaData[key].keys():
                    for mDKey, mDVal in self.dbMetaData[key][vKey].items():
                        mDVal = self.csvStr(mDVal)
                        mDValID = self.filterTree.tv.insert(vValID,
                                                            'end',
                                                            text=mDKey,
                                                            values=(mDVal,""))

    def quit(self, *args):
        """
        """
        try:
            self.savePref()
        except:
            print("savePref() error on mpdbTK exit.  Possible data loss.")
        self.master.destroy()

    def savePref(self):
        mpdb.savePref(pref=mpdb.pref)
        
    def selectItem(self, event):
        curItem = self.filterTree.tv.item(self.filterTree.tv.focus())
        col = self.filterTree.tv.identify_column(event.x)
        print ('curItem = ', curItem)
        print ('col = ', col)

        if col == '#0':
            cell_value = curItem['text']
        elif col == '#1':
            cell_value = curItem['values'][0]
        elif col == '#2':
            cell_value = curItem['values'][1]
        print ('cell_value = ', cell_value)
            
    def start(self,opts):
        """
        Call after __init__ to process command line arguments, read the
        user preferences file (or create it if it does not exist),
        and connect signals to slot functions defined below

        Parameters:
            opts, argparse dictionary.  Argparse is used in main() to create 
                  this dictionary from the command line.

        """
        self.opts = opts
        self.cmdLineArgs()
        if not self.ignoreUserPref:
            mpdb.pref = mpdb.loadPref(mpdb.prefFQPN)
            self.menubar.updateRecents(mpdb.pref['lastOpenFQPNs'],
                                       getattr(self,'openDB'))
        else:
            mpdb.pref = mpdb.defaultPref()


        if self.showSplash:
            splash = Splash(self)

        #do more init stuff here
        
        if self.showSplash:
            ## 2 second delay while loading
            sleep(1)
            ## finished loading so destroy splash
            splash.destroy()

        ## show window again
        self.master.deiconify()
        
class Splash(tk.Toplevel):
    """
    ref:
        http://code.activestate.com/recipes/577271-tkinter-splash-screen/
    """
    def __init__(self, parent, width=0.3, height=0.3, useFactor=True):
        tk.Toplevel.__init__(self, parent)
        self.title("splash")
        # get screen width and height
        ws = self.winfo_screenwidth()
        hs = self.winfo_screenheight()
        w = (useFactor and ws*width) or width
        h = (useFactor and ws*height) or height
        # calculate position x, y
        x = (ws/2) - (w/2) 
        y = (hs/2) - (h/2)
        self.geometry('%dx%d+%d+%d' % (w, h, x, y))
        
        self.overrideredirect(True)
        self.lift()
        #self.update()

def main(*args, **kwargs):
    """
    handle command line arguments and start the tk ui
    """
    pass

if __name__ == '__main__':
    """
    from ipython/jupyter qtconsole:
    %gui
    %run gui/mpdbTk.py
    """
    #ui = main()
    try:
        __IPYTHON__
        #get_ipython().magic(u'gui')
        get_ipython().magic(u'gui tk')
        root = tk.Tk()
        ui = GUI(root)
        p = argparse.ArgumentParser(description='%s: A tkinter GUI' %
                                    ui.appName)
        p.add_argument('-i', '--ignoreUserPref', action='store_true',
                   help='ignore saved user preferences',default=False)    
        p.add_argument('-s', '--showSplash', action='store_true',
                   help='display the splash screen',default=False)    
        p.add_argument('--version', action='version',
                       version='%(prog)s Version: ' + ui.__version__)
        p.add_argument('-t', '--testVal', action='store', default='',
                   metavar='cmd line arg test value',
                   help='used for development and testing')

        args = p.parse_args()

        ui.start(vars(args))

    except NameError:
        root.mainloop()
