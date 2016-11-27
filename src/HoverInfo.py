# -*- coding: utf-8 -*-
#==============================================================================
# HoverInfo.py
#------------------------------------------------------------------------------
# HoverInfo displays balloons when a text label is hovered. The original code
# can be found in the link below:
# https://jakirkpatrick.wordpress.com/2012/02/01/making-a-hovering-box-in-tkinter/
#==============================================================================

import re, sys

if sys.version_info[0] < 3:
    from Tkinter import *
else:
    from tkinter import *
 
class HoverInfo(Menu):
    
    def __init__(self, parent, text, command = None):
        self._com = command
        Menu.__init__(self, parent, tearoff = 0)
        if not isinstance(text, str):
            raise TypeError("Trying to initialise a Hover Menu with a non string type: " + text.__class__.__name__)
        toktext = re.split("\n", text)
        for t in toktext:
            self.add_command(label = t)
        self._displayed = False
        self.master.bind("<Enter>", self.Display)
        self.master.bind("<Leave>", self.Remove)
     
#    def __del__(self):
#        self.master.unbind("<Enter>")
#        self.master.unbind("<Leave>")
     
    def Display(self, event):
        if not self._displayed:
            self._displayed = True
            self.post(event.x_root, event.y_root)
        if self._com != None:
            self.master.unbind_all("<Return>")
            self.master.bind_all("<Return>", self.Click)
     
    def Remove(self, event):
        if self._displayed:
            self._displayed = False
            self.unpost()
        if self._com != None:
            self.unbind_all("<Return>")
     
    def Click(self, event):
        self._com()