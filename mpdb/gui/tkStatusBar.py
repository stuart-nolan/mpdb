"""
mpdb: Material Property Data Base as python module
tkStatusBar.py: class to define a tkinter statusbar

Revision Date: 2021.08.15

SPDX-License-Identifier: BSD-2-Clause
Copyright (c) 2021 Stuart Nolan. All rights reserved.
"""
import tkinter as tk
from tkinter import ttk
from collections import deque

class StatusBar(tk.Frame):
    """
    """
    def __init__(self, master):
        tk.Frame.__init__(self, master)
        self.label = tk.Label(self, bd=2, relief=tk.SUNKEN, anchor=tk.W)
        self.label.pack(fill=tk.X)
        self.queue = deque([])
        self.queueEmpty = True

        
    def set(self, format, *args, time=3):
        """
        Parameters:
            format, string format like "%s"
            args, strings for use in format.  i.e. "format % args"
            time, in seconds to display "format % args" in status bar
        """
        self.queue.append((format % args, time))
        if self.queueEmpty == True:
            self.display()
        
    def display(self):
        if len(self.queue) > 0:
            self.queueEmpty = False
            displayText,displayTime = self.queue.popleft()
            self.label.config(text=displayText)
            self.label.update_idletasks()
            self.label.after(displayTime*1000,self.display)
        else:
            self.queueEmpty=True
            self.clear()

    def clear(self):
        self.label.config(text="")
        self.label.update_idletasks()

