# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 14:15:00 2016

@author: William
"""
from ttk import Frame, Label, Entry
from Tkinter import *

def modelpictureframes(self,ModelPictureFrame):
        #set up the model display frame

        self.modelcanvas = Canvas(self.ModelPictureFrame,width=404,height=222)
        self.modelcanvas.pack(anchor=NW,padx=2,pady=2)
        self.modelgif = self.modelimages[0]       
        self.modelcanvas.image =self.modelgif
        #define the self call for the image currently on the canvas (to update via swapimage())
        self.image_on_canvas = self.modelcanvas.create_image(0, 0,anchor=NW, image=self.modelgif)