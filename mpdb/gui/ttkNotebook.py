"""
mpdb: Material Property Data Base as python module
ttkNotebook.py: class to define a tkinter statusbar

Revision Date: 2021.08.15

SPDX-License-Identifier: BSD-2-Clause
Copyright (c) 2021 Stuart Nolan. All rights reserved.
"""
import tkinter as tk
from tkinter import ttk

class Notebook(ttk.Notebook):
    """

    """
    def __init__(self, master, tabDefs=[{"text":"Tab0"}]):
        ttk.Notebook.__init__(self, master)
        self.tabFrames = {}
        self.tfids = {}
        for tidx,tab in enumerate(tabDefs):
            tfid = "F_%s" % tidx
            self.new_tab(tfid,tab)
        self.update_tids()
        #self.guiStyle = ttk.Style()
        #print(self.nb.tfids)
        #print(self.nb.tabFrames['F_0'].winfo_class())
        #self.guiStyle.configure('F_0.TFrame', background='wheat')
        #self.nb.tabFrames['F_0']['style']='F_0.TFrame'
        #self.guiStyle.configure('F_1.TFrame', background='#334353')
        #self.nb.tabFrames['F_1']['style']='F_1.TFrame'
        #integerTabIndex=self.index(self.tids["Tab0"])
        #pos can be one of 0, 1, ..., to maxNumTabs, or 'end'
        #pos='end'
        #self.insert(pos,ui.nb.tids["Tab0"])
        
    def get_tids(self):
        return {self.tab(v)["text"]:v for v in self.tabs()}

    def new_tab(self,tfid,tab={"text":"Tab0"}):
        self.update_tids()
        if tab["text"] not in self.tids.keys():       
            self.tabFrames.update({tfid:ttk.Frame(self)})
            self.add(self.tabFrames[tfid],**tab)
            
        if tab["text"] not in self.tfids.keys():
            self.tfids.update({tab["text"]:tfid})
            
    def rename_tab(self,tabName,newTabName):
        self.update_tids()
        
        if tabName in self.tids.keys() and newTabName not in self.tids.keys():
            self.tab(self.tids[tabName],text=newTabName)
            self.tids[newTabName] = self.tids.pop(tabName)

        if tabName in self.tfids.keys():
            self.tfids[newTabName] = self.tfids.pop(tabName)

    def update_tids(self):
        self.tids = self.get_tids()

