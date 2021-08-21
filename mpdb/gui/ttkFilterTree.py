"""
mpdb: Material Property Data Base as python module
ttkFilterTree.py: class to define a tkinter treeview that can be searched

Revision Date: 2021.08.15

SPDX-License-Identifier: BSD-2-Clause
Copyright (c) 2021 Stuart Nolan. All rights reserved.
"""
import tkinter as tk
from tkinter import ttk

class FilterTree(ttk.Frame):
    """
    """
    def __init__(self, *args, **kwargs):
        ttk.Frame.__init__(self, *args, **kwargs)
        
        self.filterKey = tk.StringVar()
        self.entry = tk.Entry(self, textvariable=self.filterKey)
        self.entry.grid(row=0, padx=4, pady=2, sticky=tk.EW)
        self.filterKey.trace("w", self.filter)
        self.detached = []
        
        self.tv = ttk.Treeview(self)
        self.tv.grid(row=1, padx=4, pady=2, sticky=tk.NSEW)
        # column and row configure needed to expand treeview in frame
        self.columnconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)

        self.tvsb = tk.Scrollbar(self)
        self.tvsb.grid(row=1,column=1, sticky=tk.NS)
        self.tv.configure(yscrollcommand=self.tvsb.set)
        self.tvsb.configure(command=self.tv.yview)
        
    def tvSortColumn(self, tvColID, reverse):
        """
        ref: https://stackoverflow.com/questions/1966929/tk-treeview-column-sort/1967793
        """
        if tvColID == '#0':
            l = [(self.tv.item(k)['text'], k) for k in \
                 self.tv.get_children('')]
        else:
            l = [(self.tv.set(k, tvColID), k) for k in \
                 self.tv.get_children('')]
        l.sort(reverse=reverse)

        # rearrange items in sorted positions
        for index, (val, k) in enumerate(l):
            self.tv.move(k, '', index)

            # reverse sort next time
            self.tv.heading(tvColID, command=lambda colID=tvColID:
                                  self.tvSortColumn(colID, not reverse))
            
    def filter(self, *args):
        pattern = self.filterKey.get()
        #avoid search on short strings
        if len(pattern) > 1:
            children = list(self.tv.get_children())
            children.extend(self.detached)
            end = len(children)
            for child in children:
                items = [self.tv.item(child, 'text')]
                items.extend([k for k in self.tv.item(child, 'values')])
                if not self.itemSearch(items,pattern):
                    if child not in self.detached:
                        self.detached.extend([child])
                    self.tv.detach(child)
                elif child in self.detached:
                    self.tv.move(child,"",end)
                    self.detached.remove(child)
                    
        elif len(pattern) < 2 and self.detached:
            children = list(self.tv.get_children())
            children.extend(self.detached)
            end = len(children)
            for item in self.detached:
                self.tv.move(item,"",end)
            self.detached = []

    def itemSearch(self, items, pattern):
        for text in items: 
            if text.lower().find(pattern.lower()) != -1:
                return True
        return False
        
