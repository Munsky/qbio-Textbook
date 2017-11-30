# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 14:06:24 2016

@author: William
"""
from ttk import Frame, Label, Entry
from Tkinter import *
def FittingsFrames(self,FittingFrame,FittingFrame2):

        FittingFrameLabel = Label(FittingFrame,text="Load Data")
        FittingFrameLabel.grid(row = 0,column = 0,padx=2,pady=2,sticky=W)
        self.OpenFile = Button(FittingFrame,text='Open Data File',command=self.openfile)
        self.OpenFile.grid(row=0,column=1,padx=2,pady=2,sticky=E)
        self.Plot_Data = Button(FittingFrame,text="Plot Data",command = self.Plot_Raw_Data)
        self.Plot_Data.grid(row=0,column=2,padx=2,pady=2)  
        
        self.Plot_Data_Hist = Button(FittingFrame, text="Plot Hist", command = self.Plot_Raw_Data_Hist)
        self.Plot_Data_Hist.grid(row=0,column=3,padx=2,pady=4)  
        
        self.FileEntry = Text(FittingFrame,width = 38,height = 2,bg=('#DCDCDC'))
        self.FileEntry.grid(row=1,column=0,padx=3,pady=2,columnspan=5)
        self.FileEntry.config(state=DISABLED)
        
        self.FittingODE = Radiobutton(FittingFrame2,text="ODE", variable=self.odeorfspvar,value=0,command=self.activate_fitting_boxes)
        self.FittingODE.grid(row=0,column=0,padx=2,pady=5)
        
        self.Fitting_FSP = Radiobutton(FittingFrame2,text="FSP", variable=self.odeorfspvar,value=1,command=self.activate_fitting_boxes)
        self.Fitting_FSP.grid(row=1,column=0,padx=2,pady=7,sticky=W)
        
        self.FittingFSP_entry = Entry(FittingFrame2,width=5,justify=RIGHT,state=NORMAL)
        self.FittingFSP_entry.grid(row=1,column=1,padx=4,pady=6)
        self.FittingFSP_entry.insert(0,"ε")
        self.FittingFSP_entry.configure(state=DISABLED)
        
        
        
        