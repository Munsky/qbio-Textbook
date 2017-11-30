# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 14:01:49 2016

@author: William
"""
from ttk import Frame, Label, Entry
from Tkinter import *
def PredictionFrames(self,PredictionsFrame):
        Predictionslabel = Label(PredictionsFrame, text="Prediction", width=7)
        Predictionslabel.grid(row=0,column=0,padx=4,pady=3,sticky=W)  
        
        self.predict = IntVar()
        self.predict.set(0)
        pre_ODEselect = Radiobutton(PredictionsFrame,text="ODE", variable=self.predict, value=0,command=self.activate_prediction_boxes)
        pre_ODEselect.grid(row=1,column=0,padx=2,pady=1,sticky=W)

        pre_SSAselect = Radiobutton(PredictionsFrame,text="SSA", variable=self.predict, value=1,command=self.activate_prediction_boxes)
        pre_SSAselect.grid(row=1, column=1,padx=2,pady=1,sticky=W)

        pre_FSPselect = Radiobutton(PredictionsFrame,text="FSP", variable=self.predict, value=2,command=self.activate_prediction_boxes)
        pre_FSPselect.grid(row=1,column = 2,padx=2,pady=1,sticky=W)
        
        self.SSA_entry = Entry(PredictionsFrame,width=7,justify=RIGHT,state=NORMAL)
        self.SSA_entry.grid(row=2,column=1,padx=5,pady=3,sticky=W)
        
        self.FSP_entry = Entry(PredictionsFrame,width=7,justify=RIGHT,state=NORMAL)
        self.FSP_entry.grid(row=2,column=2,padx=5,pady=4,sticky=W)
        
        self.SSA_entry.insert(0,"20")   
        self.FSP_entry.insert(0,"ε")
        self.FSP_entry.configure(state=DISABLED)
        self.SSA_entry.configure(state=DISABLED)
        
        
        Run_Prediction = Button(PredictionsFrame,text="Run",command=self.run_prediction)
        Run_Prediction.grid(row=2,column = 0,padx=2,pady=4)