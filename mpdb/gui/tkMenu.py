#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mpdb: Material Property Data Base as python module
tkMenu.py: class to define tkinter menus

Revision Date: 2021.08.15

SPDX-License-Identifier: BSD-2-Clause
Copyright (c) 2021 Stuart Nolan. All rights reserved.
"""

import tkinter as tk
from tkinter import ttk
import os

class tkMenubar(tk.Menu):
    """
    Usage:
        menubar = tkMenubar(parent)
        master.config(menu=menubar)
        add menu item commands via:
            menubar.menus['File'].entryconfig('Open',command=self.open)
    """
    def __init__(self, parent, *args, **kwargs):
        """
        """
        tk.Menu.__init__(self, parent, *args, **kwargs)
        self.menus = {}
        self.tearOff={"tearoff":False}
        self.makeFileMenu()
        self.makeEditMenu()
        self.makeToolsMenu()
        self.makeHelpMenu()
        
    def makeFTMenu(self):
        #FilterTree Menu
        menuDef={"label":"fTMenu"}
        menuName=menuDef["label"]
        self.menus.update({menuName:tk.Menu(self)})
        self.add_cascade(**menuDef, menu=self.menus[menuName])
        self.menus[menuName].config(**self.tearOff)
        
        #FilterTree Menu items
        cmd={"label":"Copy",
             "underline":0,
             "accelerator":"C-c"}
        self.menus[menuName].add_command(**cmd)
        subMenuDef={"label":"View Item",
                    "underline":0}
        subMenuName=subMenuDef["label"]
        self.menus.update({subMenuName:tk.Menu(self)})
        self.menus[menuName].add_cascade(**subMenuDef,
                                         menu=self.menus[subMenuName])
        self.menus[subMenuName].config(**self.tearOff)

    def makeFileMenu(self):
        #File Menu
        menuDef={"label":"File","underline":0}
        menuName=menuDef["label"]
        self.menus.update({menuName:tk.Menu(self)})
        self.add_cascade(**menuDef, menu=self.menus[menuName])
        self.menus[menuName].config(**self.tearOff)
        
        #File Menu items
        cmd={"label":"Open",
             "underline":0,
             "accelerator":"C-o"}
        self.menus[menuName].add_command(**cmd)
        subMenuDef={"label":"Open Recent",
                    "underline":5,
                    "accelerator":"C-r"}
        subMenuName=subMenuDef["label"]
        self.menus.update({subMenuName:tk.Menu(self)})
        self.menus[menuName].add_cascade(**subMenuDef,
                                         menu=self.menus[subMenuName])
        self.menus[subMenuName].config(**self.tearOff)
        cmd={"label":"Save",
             "underline":0,
             "accelerator":"C-s"}
        self.menus[menuName].add_command(**cmd)
        self.menus[menuName].add_separator()
        cmd={"label":"Quit",
             "underline":0,
             "accelerator":"C-q"}
        self.menus[menuName].add_command(**cmd)

    def makeEditMenu(self):
        #Edit Menu
        menuDef={"label":"Edit","underline":0}
        menuName=menuDef["label"]
        self.menus.update({menuName:tk.Menu(self)})
        self.add_cascade(**menuDef, menu=self.menus[menuName])
        self.menus[menuName].config(**self.tearOff)
        
        #Edit Menu items
        cmd={"label":"Copy",
             "underline":0,
             "accelerator":"C-c"}
        self.menus[menuName].add_command(**cmd)

    def makeToolsMenu(self):
        #Tools Menu
        menuDef={"label":"Tools","underline":0}
        menuName=menuDef["label"]
        self.menus.update({menuName:tk.Menu(self)})
        self.add_cascade(**menuDef, menu=self.menus[menuName])
        self.menus[menuName].config(**self.tearOff)

        #Edit Menu items
        cmd={"label":"Clear Open Recent"}
        self.menus[menuName].add_command(**cmd)

    def makeHelpMenu(self):
        #Help Menu
        menuDef={"label":"Help","underline":0}
        menuName=menuDef["label"]
        self.menus.update({menuName:tk.Menu(self)})
        self.add_cascade(**menuDef, menu=self.menus[menuName])
        self.menus[menuName].config(**self.tearOff)
        
    def updateRecents(self,fileList,openFunc):
        """
        update the 'Open Recent' menu
        
        Parameters:
            fileList list of FQPNs for each recent file (i.e. self.pref['
        """
        self.clearRecents()
        for idx,item in enumerate(fileList):
            displayName=str(idx)+": "+os.path.split(item)[1]
            self.menus['Open Recent'].add_command(label=displayName,
                                                  command=lambda f=item: openFunc(f))

    def clearRecents(self):
        """
        """
        lastRecent = self.menus['Open Recent'].index('end')
        self.menus['Open Recent'].delete(0,lastRecent)

