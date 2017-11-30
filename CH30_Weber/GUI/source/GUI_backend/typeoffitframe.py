# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 14:10:25 2016

@author: William

Python File that creates the fitting frame and buttons for the gui
"""
from ttk import Frame, Label, Entry
from Tkinter import *
def TypeOfFittingFrames(self,TypeofFitFrame):

        self.fittype = IntVar()
        self.fittype.set(0)        
        FminFrame = Frame(TypeofFitFrame)
        FminFrame.grid(row=0,column=0)
        TypeFitLabel = Label(FminFrame,text="Parameter Search")
        TypeFitLabel.grid(row=0,column=0,padx=2,pady=2,sticky=W)
        
        OptimizationFmin = Radiobutton(FminFrame,text="scipy.minimize    ",variable = self.fittype, value =0,command=self.activate_fittype_box)
        OptimizationFmin.grid(row=1,column=0,padx=2,pady=2,columnspan = 2)
        self.fminmethod = StringVar()
        self.Fmin_Method = OptionMenu(FminFrame,self.fminmethod,"L-BFGS-B", "    TNC    ","   SLSQP   ") 
        self.fminmethod.set('L-BFGS-B')
        self.Fmin_Method.grid(row=1,column=2,padx=2,pady=2,columnspan = 3)
        
        SAFrame = Frame(TypeofFitFrame)
        SAFrame.grid(row=1,column=0)
        
        
        OptimizationSA = Radiobutton(SAFrame,text="SA              ",variable = self.fittype, value =1,command=self.activate_fittype_box)
        OptimizationSA.grid(row=2,column=0,padx=2,pady=2)
        SA_label = Label(SAFrame,text="Max Iterations")
        SA_label.grid(row=2,column=1,columnspan=3)
        self.SA_limit = Entry(SAFrame,width=7,justify=RIGHT)
        self.SA_limit.grid(row=2,column=4,padx=3,pady=2,columnspan = 2,sticky=W)
        self.SA_limit.insert(0,"500")
        self.SA_limit.config(state=DISABLED)
        
        
        OptimizationMetHaste = Radiobutton(TypeofFitFrame,text="Met-Haste",variable = self.fittype, value =2,command=self.activate_fittype_box)
        OptimizationMetHaste.grid(row=3,column=0,padx=2,pady=2)
        
        self.burn_label = Label(TypeofFitFrame, text="Burn")
        self.burn_label.grid(row=3,column=1,padx=2,pady=2)
        
        self.chain_label = Label(TypeofFitFrame, text="Chain")
        self.chain_label.grid(row=3,column=3,padx=1,pady=2)
        
        self.thin_label = Label(TypeofFitFrame, text="Thin")
        self.thin_label.grid(row=3,column=5,padx=2,pady=2)
        
        self.thin_label = Label(TypeofFitFrame, text="Mut R")
        self.thin_label.grid(row=3,column=7,padx=2,pady=2)        
        
        self.N_Burn = Entry(TypeofFitFrame,width=5,justify=RIGHT,state=NORMAL)
        self.N_Burn.grid(row=3,column=2,padx=3,pady=2)
        
        self.N_Chain = Entry(TypeofFitFrame,width=5,justify=RIGHT,state=NORMAL)
        self.N_Chain.grid(row=3,column=4,padx=3,pady=2)
        
        self.N_Thin = Entry(TypeofFitFrame,width=5,justify=RIGHT,state=NORMAL)
        self.N_Thin.grid(row=3,column=6,padx=3,pady=2)
        
        self.Mut_Rate = Entry(TypeofFitFrame,width=5,justify=RIGHT,state=NORMAL)
        self.Mut_Rate.grid(row=3,column=8,padx=4,pady=2)        
        
        Plot_Data = Button(TypeofFitFrame,text="Plot Data",command = self.Plot_Raw_Data)
        Plot_Data.grid(row=1,column=7,padx=4,pady=2,columnspan=2)         
        
        self.N_Burn.insert(0,"10")   
        self.N_Thin.insert(0,"10")
        self.N_Chain.insert(0,"100") 
        self.Mut_Rate.insert(0,".2")
        
        self.N_Burn.config(state=DISABLED)   
        self.N_Thin.config(state=DISABLED)
        self.N_Chain.config(state=DISABLED)   
        self.Mut_Rate.config(state=DISABLED)  
        
        
        self.MetHasteParSpace = Button(TypeofFitFrame,text="View Parameter Space",state=DISABLED,command=self.PlotParameterSpace)
        self.MetHasteParSpace.grid(row=4,column=5,padx=3,pady=3,columnspan=4)
        
        checkboxframe = Frame(TypeofFitFrame)
        checkboxframe.grid(row=4,column=0,columnspan=5)
        self.k12checkvar = IntVar()
        self.k12checkvar.set(0)
        self.k12check = Checkbutton(checkboxframe,text="k12",variable=self.k12checkvar)
        self.k12check.grid(row=0,column=0,padx=2,pady=2)

        self.k23checkvar = IntVar()
        self.k23checkvar.set(0)
        self.k23check = Checkbutton(checkboxframe,text="k23",variable=self.k23checkvar)
        self.k23check.grid(row=0,column=1,padx=2,pady=2)
        
        self.k21checkvar = IntVar()
        self.k21checkvar.set(0)
        self.k21check = Checkbutton(checkboxframe,text="k21",variable=self.k21checkvar)
        self.k21check.grid(row=0,column=2,padx=2,pady=2)

        self.k32checkvar = IntVar()
        self.k32checkvar.set(0)
        self.k32check = Checkbutton(checkboxframe,text="k32",variable=self.k32checkvar)
        self.k32check.grid(row=1,column=0,padx=2,pady=2)
    
        self.kr2checkvar = IntVar()
        self.kr2checkvar.set(0)
        self.kr2check = Checkbutton(checkboxframe,text="kr2",variable=self.kr2checkvar)
        self.kr2check.grid(row=1,column=1,padx=2,pady=2)

        self.kr3checkvar = IntVar()
        self.kr3checkvar.set(0)
        self.kr3check = Checkbutton(checkboxframe,text="kr3",variable=self.kr3checkvar)
        self.kr3check.grid(row=1,column=2,padx=2,pady=2)
        
        self.gmmcheckvar = IntVar()
        self.gmmcheckvar.set(0)
        self.gmmcheck = Checkbutton(checkboxframe,text=" γ ",variable=self.gmmcheckvar)
        self.gmmcheck.grid(row=1,column=3,padx=2,pady=2)

        self.betcheckvar = IntVar()
        self.betcheckvar.set(0)
        self.betcheck = Checkbutton(checkboxframe,text=" β ",variable=self.betcheckvar)
        self.betcheck.grid(row=0,column=3,padx=2,pady=2)
        
        self.k12check.configure(state=DISABLED)
        self.k21check.configure(state=DISABLED)
        self.k23check.configure(state=DISABLED)
        self.k32check.configure(state=DISABLED)
        self.kr2check.configure(state=DISABLED)
        self.kr3check.configure(state=DISABLED)
        self.gmmcheck.configure(state=DISABLED)
        self.betcheck.configure(state=DISABLED)


        
        FitRun = Button(TypeofFitFrame,text="  Try Fit  ",command=self.Run_Fit)
        FitRun.grid(row=5,column=0,padx=2,pady=3,columnspan=1)
        
        UpdatePar = Button(TypeofFitFrame,text="Update Parameters",command=self.Update_Parameters)
        UpdatePar.grid(row=5,column=1,padx=2,pady=3,columnspan=3)
        
        SavePar = Button(TypeofFitFrame,text="Save Parameters")
        SavePar.grid(row=5,column=4,padx=2,pady=3,columnspan=3)
        
        Moments = Button(TypeofFitFrame,text="2nd Moments",command=self.Second_Moment_Call)
        Moments.grid(row=5,column=7,padx=5,pady=3,columnspan=2)
