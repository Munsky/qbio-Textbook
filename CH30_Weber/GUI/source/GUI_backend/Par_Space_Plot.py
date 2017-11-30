# -*- coding: utf-8 -*-
"""
Created on Sat Jun 18 17:01:49 2016

@author: William

Python file that does the iterations and loops for the parameter space plots 2x2 - 7x7 var
"""
from ttk import Frame, Label, Entry
from Tkinter import *
from matplotlib import cm
from numpy import linspace
def ParSpace_Plot(self,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist,**kwargs):
    
    
    a = kwargs.get('alpha',1)
    if a =="alpha":
        a=1
    color = kwargs.get('color','Solid')
    
    
    #################################################################
    # Preallocate variables based on how many checkboxes were checked
    # Preallocate lists to plot as well
    #################################################################
    if sum([k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot]) == 2:
            var1 = [] #first variable 
            var2 = []
            used = [] #used to not reiterate variables
            var1,var1str,used = self.par_space_setup(var1,used,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist)
            var2,var2str,used = self.par_space_setup(var2,used,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist)
            varlist = [var1,var2]
            varstrlist = [var1str,var2str] 
            cornerindex = 3 #index of the corner plot
            xindex = [4] #index for which plots should have labels x 
            yindex = [1] #index for which plots should have labels y

            
    if sum([k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot]) == 3:
            
            var1 = []
            var2 = []
            var3 = []
            used = []
            var1,var1str,used = self.par_space_setup(var1,used,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist)
            var2,var2str,used = self.par_space_setup(var2,used,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist)
            var3,var3str,used = self.par_space_setup(var3,used,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist)
            varlist = [var1,var2,var3]
            varstrlist = [var1str,var2str,var3str]
            cornerindex = 7
            xindex = [8,9]
            yindex = [1,4]



    if sum([k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot]) == 4:
            
            var1 = []
            var2 = []
            var3 = []
            var4 = []
            used = []
            var1,var1str,used = self.par_space_setup(var1,used,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist)
            var2,var2str,used = self.par_space_setup(var2,used,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist)
            var3,var3str,used = self.par_space_setup(var3,used,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist)
            var4,var4str,used = self.par_space_setup(var4,used,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist)
            varlist = [var1,var2,var3,var4]
            varstrlist = [var1str,var2str,var3str,var4str]
            cornerindex = 13
            xindex = [14,15,16]
            yindex = [1,5,9]


    if sum([k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot]) == 7:
            
            var1 = []
            var2 = []
            var3 = []
            var4 = []
            var5=[]
            var6=[]
            var7 = []
            used = []
            var1,var1str,used = self.par_space_setup(var1,used,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist)
            var2,var2str,used = self.par_space_setup(var2,used,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist)
            var3,var3str,used = self.par_space_setup(var3,used,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist)
            var4,var4str,used = self.par_space_setup(var4,used,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist)
            var5,var5str,used = self.par_space_setup(var5,used,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist)
            var6,var6str,used = self.par_space_setup(var6,used,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist)
            var7,var7str,used = self.par_space_setup(var7,used,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist)
            varlist = [var1,var2,var3,var4,var5,var6,var7]
            varstrlist = [var1str,var2str,var3str,var4str,var5str,var6str,var7str]
            cornerindex = 43
            xindex = [44,45,46,47,48,49]
            yindex = [1,8,15,22,29,36]

    if sum([k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot]) == 8:
            
            var1 = []
            var2 = []
            var3 = []
            var4 = []
            var5=[]
            var6=[]
            var7 = []
            var8 = []
            used = []
            var1,var1str,used = self.par_space_setup(var1,used,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist)
            var2,var2str,used = self.par_space_setup(var2,used,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist)
            var3,var3str,used = self.par_space_setup(var3,used,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist)
            var4,var4str,used = self.par_space_setup(var4,used,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist)
            var5,var5str,used = self.par_space_setup(var5,used,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist)
            var6,var6str,used = self.par_space_setup(var6,used,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist)
            var7,var7str,used = self.par_space_setup(var7,used,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist)
            var8,var8str,used = self.par_space_setup(var8,used,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist)
            
            varlist = [var1,var2,var3,var4,var5,var6,var7,var8]
            varstrlist = [var1str,var2str,var3str,var4str,var5str,var6str,var7str,var8str]
            cornerindex = 57
            xindex = [58,59,60,61,62,63,64]
            yindex = [1,9,17,25,33,41,49]
    
    if sum([k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot]) == 6:
            
            var1 = []
            var2 = []
            var3 = []
            var4 = []
            var5=[]
            var6=[]
            used = []
            var1,var1str,used = self.par_space_setup(var1,used,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist)
            var2,var2str,used = self.par_space_setup(var2,used,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist)
            var3,var3str,used = self.par_space_setup(var3,used,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist)
            var4,var4str,used = self.par_space_setup(var4,used,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist)
            var5,var5str,used = self.par_space_setup(var5,used,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist)
            var6,var6str,used = self.par_space_setup(var6,used,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist)
            varlist = [var1,var2,var3,var4,var5,var6]
            varstrlist = [var1str,var2str,var3str,var4str,var5str,var6str]
            cornerindex = 31
            xindex = [32,33,34,35,36]
            yindex = [1,7,13,19,25]

        
    if sum([k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot]) == 5:
            
            var1 = []
            var2 = []
            var3 = []
            var4 = []
            var5=[]
            used = []
            var1,var1str,used = self.par_space_setup(var1,used,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist)
            var2,var2str,used = self.par_space_setup(var2,used,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist)
            var3,var3str,used = self.par_space_setup(var3,used,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist)
            var4,var4str,used = self.par_space_setup(var4,used,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist)
            var5,var5str,used = self.par_space_setup(var5,used,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist)
            varlist = [var1,var2,var3,var4,var5]
            varstrlist = [var1str,var2str,var3str,var4str,var5str]
            cornerindex = 21
            xindex = [22,23,24,25]
            yindex = [1,6,11,16]
            
                                                                        
            
    k = 0

    time = linspace(0,1,len(varlist[0][:]))
        
    #plot the variables vs each other 
    for i in range(0,len(varlist)):
        for j in range(0,len(varlist)):
            k+=1
        
            if i == j: # if the variable is being compared against itself, use a histogram
                subplot = self.fig.add_subplot(len(varlist),len(varlist),k)
                
                if color=="Greys":
                    
                    subplot.hist(varlist[i],facecolor="grey")  
                else:
                    subplot.hist(varlist[i])
                
            if i > j: #if its below the diagonal plot
            
                subplot = self.fig.add_subplot(len(varlist),len(varlist),k)
                
                
                if color == "Jet":                
                    subplot.scatter(varlist[j],varlist[i],c= time,cmap='jet',marker='.',alpha=a,linewidths=None,edgecolors='face')
                    
                if color == "Solid":
                    subplot.scatter(varlist[j],varlist[i],marker='.',alpha=a,linewidths=None,edgecolors='face')

                if color == "Greys":
                    subplot.scatter(varlist[j],varlist[i],c= time,cmap='Greys',marker='.',alpha=a,linewidths=None,edgecolors='face')
                    
                if color == "Tri":
                    subplot.scatter(varlist[j][0:len(varlist[j])/3],varlist[i][0:len(varlist[i])/3],c='r',marker='.',alpha=a,linewidths=None,edgecolors='face')
                    subplot.scatter(varlist[j][len(varlist[j])/3:2*len(varlist[j])/3],varlist[i][len(varlist[i])/3:2*len(varlist[j])/3],c='b',marker='.',alpha=a,linewidths=None,edgecolors='face')
                    subplot.scatter(varlist[j][2*len(varlist[j])/3:len(varlist[j])],varlist[i][2*len(varlist[i])/3:len(varlist[j])],c='y',marker='.',alpha=a,linewidths=None,edgecolors='face')
                    
                subplot.get_xaxis().set_visible(False)
                subplot.get_yaxis().set_visible(False)
                
           
            #plot the corner plots with both x and y labels
            if k == cornerindex-1:
                subplot.get_xaxis().set_visible(True)
                subplot.set_xlabel('$'+varstrlist[j] +'$')
            if k == cornerindex:
                subplot.get_xaxis().set_visible(True)
                subplot.get_yaxis().set_visible(True)
                subplot.set_xlabel('$'+varstrlist[j]+'$')
                subplot.set_ylabel('$'+varstrlist[i]+'$')
            #plot anything on left border with y labels
            if k in yindex:
                subplot.get_yaxis().set_visible(True)
                subplot.set_ylabel('$'+varstrlist[i]+'$')
            #plot anything on the bottom border with x labels
            if k in xindex:
                subplot.get_xaxis().set_visible(True)
                subplot.set_xlabel('$'+varstrlist[j]+'$')
            if i == j:
                subplot.get_yaxis().set_visible(True)
                
    self.figurecanvas.show()
    self.fig.canvas.draw()

 