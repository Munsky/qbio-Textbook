# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 13:58:08 2016

@author: William
"""
from ttk import Frame, Label, Entry
from Tkinter import *
def x0Frames(self,x0Frame):
    
        x0label = Label(self.x0Frame,text="Initial Values",width = 10)
        x0label.pack(anchor=NW,padx=2,pady=3)
        
        S1_label = Label(self.x0Frame,text="S1",width=3)  
        S1_label.pack(side=LEFT,pady=3)
        
        self.x0_S1 = Entry(self.x0Frame,text="S1",width=6,justify=RIGHT)
        self.x0_S1.pack(side = LEFT,padx=1,pady=3)
        self.x0_S1.insert(0,"1")
        
        S2_label = Label(self.x0Frame,text="S2",width=3)  
        S2_label.pack(side=LEFT,pady=3)
        
        self.x0_S2 = Entry(self.x0Frame,text="S2",width=6,justify=RIGHT)
        self.x0_S2.pack(side = LEFT,padx=1,pady=3)
        self.x0_S2.insert(0,"0")

        S3_label = Label(self.x0Frame,text="S3",width=3)  
        S3_label.pack(side=LEFT,padx=1,pady=3)
        
        self.x0_S3 = Entry(self.x0Frame,text="S3",width=6,justify=RIGHT)
        self.x0_S3.pack(side = LEFT, padx=1,pady=3)
        self.x0_S3.insert(0,"0")
        
        R_label = Label(self.x0Frame,text="R",width=2)  
        R_label.pack(side=LEFT,padx=1,pady=3)

        self.x0_R = Entry(self.x0Frame,text="R",width=6,justify=RIGHT)
        self.x0_R.pack(side=LEFT,padx=1,pady=3)
        self.x0_R.insert(0,"0")
        
        t0_label = Label(self.x0Frame,text="t0",width=2)  
        t0_label.pack(side=LEFT,padx=1,pady=3)

        self.x0_t0 = Entry(self.x0Frame,text="t0",width=6,justify=RIGHT)
        self.x0_t0.pack(side=LEFT,padx=1,pady=3)
        self.x0_t0.insert(0,"0")

        tf_label = Label(self.x0Frame,text="tf",width=2)  
        tf_label.pack(side=LEFT,padx=1,pady=3)

        self.x0_tf = Entry(self.x0Frame,text="tf",width=6,justify=RIGHT)
        self.x0_tf.pack(side=LEFT,padx=1,pady=3)
        self.x0_tf.insert(0,"100")
        
        space_label2 = Label(self.x0Frame,text="",width = 0)
        space_label2.pack(side=LEFT,padx=3,pady=3)
#_____________________________________________________________________________