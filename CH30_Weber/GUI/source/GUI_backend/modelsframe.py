# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 11:34:19 2016

@author: William
"""

from ttk import Frame, Label, Entry
from Tkinter import *
def ModelsInputsFrames(self,ModelsFrame,InputFrame):
        #set up the model selection frame

        
        #define all the radio buttons and label them
        methodlabel = Label(ModelsFrame, text="Model:", width=6)
        methodlabel.pack(anchor=NW, padx=2, pady=2)  
        
        modelselect0 = Radiobutton(ModelsFrame,text="0", variable=self.case, value=0, command = self.swapimage)
        modelselect0.pack(side=LEFT)
        
        modelselect1 = Radiobutton(ModelsFrame,text="1", variable=self.case, value=1,command = self.swapimage)
        modelselect1.pack(side=LEFT)
        
        modelselect2 = Radiobutton(ModelsFrame,text="2", variable=self.case, value=2,command = self.swapimage)
        modelselect2.pack(side=LEFT)

        modelselect3 = Radiobutton(ModelsFrame,text="3", variable=self.case, value=3,command = self.swapimage)
        modelselect3.pack(side=LEFT)
        
        modelselect4 = Radiobutton(ModelsFrame,text="4", variable=self.case, value=4,command = self.swapimage)
        modelselect4.pack(side=LEFT,padx=2,pady=2)

        #define the input frame
      
        
        Inputlabel = Label(InputFrame, text="Input F(t)                                                          ")
        Inputlabel.pack(anchor=NW,padx=4,pady=2)      
        
        self.Input_entry = Entry(InputFrame,width=37,justify=RIGHT)
        self.Input_entry.pack(side=RIGHT,padx=4,pady=5)
        self.Input_entry.insert(0,"(1-cos(2*3.14/30*t))*(t>5)*(t<70)")    
        