# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 13:52:15 2016

@author: William

python file that creates the variable entries for the GUI
"""
from ttk import Frame, Label, Entry
from Tkinter import *
def VariablesFrames(self,variablesFrame):

        
        self.K12Frame= Frame(self.variablesFrame)
        self.K12Frame.grid(row=1,column=0,padx = 2,pady=2,rowspan=1)        

        self.varFrame= Frame(self.variablesFrame)
        self.varFrame.grid(row=0,column=0,padx = 1,pady=1,rowspan=1)    
        
        self.variablelabel = Label(self.varFrame, text="Parameters", width=10)
        self.variablelabel.pack(anchor=NW)  
        
        self.k12_label = Label(self.K12Frame, text="k" +u"\u2081"+u"\u2082"+" ")
        self.k12_label.pack(side=LEFT,padx=1)
        self.k12_entry = Entry(self.K12Frame,width=7,justify=RIGHT)
        self.k12_entry.pack(side=RIGHT,padx=2,pady=2)
        self.k12_entry.insert(0,"1")        
        
        self.K23Frame= Frame(self.variablesFrame)
        self.K23Frame.grid(row=2,column=0,padx = 2,pady=2)              
             
        
        self.k23_label = Label(self.K23Frame, text="k" +u"\u2082"+u"\u2083"+" ")
        self.k23_label.pack(side=LEFT)
        self.k23_entry = Entry(self.K23Frame,width=7,justify=RIGHT)
        self.k23_entry.pack(side=RIGHT,padx=2,pady=3)
        self.k23_entry.insert(0,"1")  

        self.K21Frame= Frame(self.variablesFrame)
        self.K21Frame.grid(row=0,column=1,padx = 2,pady=2) 
        
        self.k21_label = Label(self.K21Frame, text="k"+ u"\u2082"+ u"\u2081"+" ")
        self.k21_label.pack(side=LEFT)
        self.k21_entry = Entry(self.K21Frame,width=7,justify=RIGHT)
        self.k21_entry.pack(side=RIGHT,padx=2,pady=3)
        self.k21_entry.insert(0,".5")  
        
        self.K32Frame= Frame(self.variablesFrame)
        self.K32Frame.grid(row=1,column=1,padx = 2,pady=2) 
        
        self.k32_label = Label(self.K32Frame, text="k"+ u"\u2083"+ u"\u2082"+" ")
        self.k32_label.pack(side=LEFT)
        self.k32_entry = Entry(self.K32Frame,width=7,justify=RIGHT)
        self.k32_entry.pack(side=RIGHT,padx=2,pady=2)
        self.k32_entry.insert(0,".5")  
        
        self.Kr2Frame= Frame(self.variablesFrame)
        self.Kr2Frame.grid(row=2,column=1,padx = 2,pady=2) 
        
        self.kr2_label = Label(self.Kr2Frame, text="k"+ u"\u1D63"+ u"\u2082"+" ")
        self.kr2_label.pack(side=LEFT)
        self.kr2_entry = Entry(self.Kr2Frame,width=7,justify=RIGHT)
        self.kr2_entry.pack(side=RIGHT,padx=3,pady=4)
        self.kr2_entry.insert(0,"2")  
        
        self.Kr3Frame= Frame(self.variablesFrame)
        self.Kr3Frame.grid(row=0,column=2,padx = 2,pady=2) 
        
        self.kr3_label = Label(self.Kr3Frame, text="k"+ u"\u1D63"+ u"\u2083"+" ")
        self.kr3_label.pack(side=LEFT)
        self.kr3_entry = Entry(self.Kr3Frame,width=7,justify=RIGHT)
        self.kr3_entry.pack(side=RIGHT,padx=3,pady=2)
        self.kr3_entry.insert(0,".5")        
        
        self.gammaFrame= Frame(self.variablesFrame)
        self.gammaFrame.grid(row=1,column=2,padx = 2,pady=2) 
        
        self.gamma_label = Label(self.gammaFrame, text=" γ  ")
        self.gamma_label.pack(side=LEFT)
        self.gamma_entry = Entry(self.gammaFrame,width=7,justify=RIGHT)
        self.gamma_entry.pack(side=RIGHT,padx=3,pady=2)
        self.gamma_entry.insert(0,".2")
        
        self.betaFrame= Frame(self.variablesFrame)
        self.betaFrame.grid(row=2,column=2,padx = 2,pady=2) 
        
        self.beta_label = Label(self.betaFrame, text=" β ")
        self.beta_label.pack(side=LEFT)
        self.beta_entry = Entry(self.betaFrame,width=7,justify=RIGHT)
        self.beta_entry.pack(side=RIGHT,padx=3,pady=2)
        self.beta_entry.insert(0,"1")        
        