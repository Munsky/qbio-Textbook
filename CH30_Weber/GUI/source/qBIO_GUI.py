# -*- coding: utf-8 -*-
"""
Created on Tue May 24 10:58:11 2016

@author: William Raymond
wsraymon@rams.colostate.edu

"""


##############################################################################
# CURRENT BUGS AND PROBLEMS

# no FSP fit
# no about section

# placeholders for menubar functions
##############################################################################

from Tkinter import * #import tkinter 
import os
from tkColorChooser import askcolor  

import tkFont
from ttk import Frame, Label,Entry, Progressbar, Style #import necessary tkk functions


from tkFileDialog import askopenfilename, asksaveasfilename #tkinter file opening dialouge


import numpy as np #numpy
import tkMessageBox 
from scipy.integrate import odeint 
import tkFileDialog #open/save dialog
import time
import pickle

from matplotlib.figure import Figure #figure handler
from matplotlib import cm
import matplotlib.pyplot as plt #plot

from scipy.integrate import ode

from scipy import linalg
from scipy.optimize import minimize
from scipy.interpolate import spline    
from scipy import sparse





class GeneModel():    
    
    def __init__(self):
        pass

    def Create_CME_Sparse(self,parameters,N):
    
        k12 = parameters[0] #define parameters
        k23 = parameters[1]
        k21 = parameters[2]
        k32 = parameters[3]
        kr2 = parameters[4]
        kr3 = parameters[5]
        gamma = parameters[6]

        A = sparse.csr_matrix((3*N,3*N))  #call scipy's sparse to create the matrix
        
        maindiag = np.array([[-k12, (-k23-k21-kr2), (-k32-kr3)]])  #make the middle diagonal
        
        maindiag = np.tile(maindiag,N)
        
        zerosdiag = np.array([[0,0,0]])
        zerosdiag = np.tile(zerosdiag,N) #any diagonal that is not used
        
        downoffdiag = np.array([[k12, k23,0]])  #down one diag 
        downoffdiag = np.tile(downoffdiag,N)  #tile = repeat the vector N times

        
        upoffdiag = np.array([[k21, k32,0]])  #Up one diagonal
        upoffdiag = np.tile(upoffdiag,N)
        
        degradediag = np.array([[gamma, gamma, gamma]])
        maindegradediag = np.array([[0,0,0]])
        for i in range(1,N):
            degradediag = np.append(degradediag,[(i+1)*gamma, (i+1)*gamma,(i+1)*gamma])
            maindegradediag = np.append(maindegradediag,[(i)*gamma,(i)*gamma,(i)*gamma])  #degrade diagonals
            
        degradediag = np.array([degradediag])      
        #maindegradediag[0,0:3] = maindegradediag[0,0:3]+1        
        maindiag = maindiag-maindegradediag

        #maindiag[0][0:3] = maindiag[0][0:3] + [0, kr2, kr3]
    
    
            
        transcriptdiag = np.array([[0, kr2,kr3]])
        transcriptdiag = np.tile(transcriptdiag,N-1)   #transcriptions 
        
        A = np.zeros((3*N,3*N))
        np.fill_diagonal(A,maindiag)
        np.fill_diagonal(A[:-1, 1:], upoffdiag)
        np.fill_diagonal(A[1:, :-1], downoffdiag)
        np.fill_diagonal(A[:-3, 3:], degradediag)
        np.fill_diagonal(A[3:,:-3], transcriptdiag)  #place into matrix 
        return A
        

    
    def Create_CME(self,parameters,N):
    
        k12 = parameters[0]
        k23 = parameters[1]
        k21 = parameters[2]
        k32 = parameters[3]
        kr2 = parameters[4]
        kr3 = parameters[5]
        
        gamma = parameters[6]
        A = np.zeros((3*N,3*N))
        
        maindiag = np.array([[-k12, (-k23-k21-kr2), (-k32-kr3)]])
        
        maindiag = np.tile(maindiag,N)
        
        zerosdiag = np.array([[0,0,0]])
        zerosdiag = np.tile(zerosdiag,N)
        
        downoffdiag = np.array([[k12, k23,0]])
        downoffdiag = np.tile(downoffdiag,N)

        
        upoffdiag = np.array([[k21, k32,0]])
        upoffdiag = np.tile(upoffdiag,N)
        
        degradediag = np.array([[gamma, gamma, gamma]])
        maindegradediag = np.array([[0,0,0]])
        for i in range(1,N):
            degradediag = np.append(degradediag,[(i+1)*gamma, (i+1)*gamma,(i+1)*gamma])
            maindegradediag = np.append(maindegradediag,[(i)*gamma,(i)*gamma,(i)*gamma])
            
        degradediag = np.array([degradediag])      
        #maindegradediag[0,0:3] = maindegradediag[0,0:3]+1        
        maindiag = maindiag-maindegradediag

        #maindiag[0][0:3] = maindiag[0][0:3] + [0, kr2, kr3]
    
    
            
        transcriptdiag = np.array([[0, kr2,kr3]])
        transcriptdiag = np.tile(transcriptdiag,N-1)
        
        A = np.zeros((3*N,3*N))
        np.fill_diagonal(A,maindiag)
        np.fill_diagonal(A[:-1, 1:], upoffdiag)
        np.fill_diagonal(A[1:, :-1], downoffdiag)
        np.fill_diagonal(A[:-3, 3:], degradediag)
        np.fill_diagonal(A[3:,:-3], transcriptdiag)
        return A

    def Create_CME_TimeVar(self,parameters,N,TimeVar,case,t):
        
        k12 = parameters[0]
        k23 = parameters[1]
        k21 = parameters[2]
        k32 = parameters[3]
        kr2 = parameters[4]
        kr3 = parameters[5]
        gamma = parameters[6]
        beta = parameters[7]        
        
        if case == 1:
            k12 = max(parameters[0] +beta*TimeVar(t),0)
        if case == 2:
            k23 = max(0,parameters[1] + beta*TimeVar(t))
        if case == 3:
            k21 = max(0,parameters[2]+ beta*TimeVar(t))
        if case == 4:
            k32 = max(0,parameters[3] + beta*TimeVar(t))
        
       
        
        H = np.array([[-k12 , k21 , 0 ],
                    [k12,-k21-k23-kr2, k32],
                    [0,    k23, -k32-kr3]])
        A = np.zeros((3*N,3*N))
        
        
        maindiag = np.array([[-k12, (-k23-k21-kr2), (-k32-kr3)]])
        
        maindiag = np.tile(maindiag,N)
        
        zerosdiag = np.array([[0,0,0]])
        zerosdiag = np.tile(zerosdiag,N)
        
        downoffdiag = np.array([[k12, k23,0]])
        downoffdiag = np.tile(downoffdiag,N)

        
        upoffdiag = np.array([[k21, k32,0]])
        upoffdiag = np.tile(upoffdiag,N)
        
        degradediag = np.array([[gamma, gamma, gamma]])
        maindegradediag = np.array([[0,0,0]])
        for i in range(1,N):
            degradediag = np.append(degradediag,[(i+1)*gamma, (i+1)*gamma,(i+1)*gamma])
            maindegradediag = np.append(maindegradediag,[(i)*gamma,(i)*gamma,(i)*gamma])
            
        degradediag = np.array([degradediag])
        
        #maindegradediag[0,0:3] = maindegradediag[0,0:3]+1        
        maindiag = maindiag-maindegradediag

        
        #maindiag[0][0:3] = maindiag[0][0:3] + [0, kr2, kr3]
    
    
        
        transcriptdiag = np.array([[0, kr2,kr3]])
        transcriptdiag = np.tile(transcriptdiag,N-1)
        
        A = np.zeros((3*N,3*N))
        np.fill_diagonal(A,maindiag)
        np.fill_diagonal(A[:-1, 1:], upoffdiag)
        np.fill_diagonal(A[1:, :-1], downoffdiag)
        np.fill_diagonal(A[:-3, 3:], degradediag)
        np.fill_diagonal(A[3:,:-3], transcriptdiag)
        return A
    
    
    
        
    def FSPCalc(self,t0,tf,A,x0,N,tol,parameters):
        ################################################################            
        #Computes the Finite State Projection for time independent systems
        #within given error.
        #A is the state reaction matrix, tf is the final time of the projections, 
        #X_0 is the initial states, error is the amount of error between the FSP and 
        #the actual solution.    
        ################################################################
        
        ones = np.array([(np.ones(np.shape(A[1])))])   #computes a 2D row of 1's.
        matexp = linalg.expm(A*tf)    
                #compute matrix exponential
        
        G = np.dot((np.dot(ones,matexp)),x0) 
        
        if G >= 1-tol:                 
            xf= np.dot(matexp,x0)       
            return xf
        else: 
            N += 10
            if N < 300:
                
                x0 =np.zeros((1,3*N)).T
                x0[0][0] = 1
                
                A = self.Create_CME(parameters,N)
                
                
                self.FSPCalc(t0,tf,A,x0,N,tol,parameters)
            else:
                return
   
            
    #Used, makes the A(t) for each given time point for the integrator to integrate from A(t1) to A(t2)
    
    def Test_A(self,parameters,N,TimeVar,case,t):
        A0 = self.Create_CME_Sparse(parameters,N)
        At = np.empty((1,len(t)),dtype=object)
        At[0][0] = A0
        
        newpar = np.copy(parameters)
        for i in range(1,len(t)):
            if case == 1:
                k12 = max(parameters[0] +parameters[7]*TimeVar(t[i]),0)
                newpar[0] = k12
            if case == 2:
                k23 = max(0,parameters[1] + parameters[7]*TimeVar(t[i]))
                newpar[1] = k23
            if case == 3:
                k21 = max(0,parameters[2]+ parameters[7]*TimeVar(t[i]))
                newpar[2] = k21
            if case == 4:
                k32 = max(0,parameters[3] + parameters[7]*TimeVar(t[i])) 
                newpar[3] = k32

            At[0][i] = self.Create_CME_Sparse(newpar,N)
            
        
        return At
                    
   
    #integrate the preallocated FSP
    def FSPCalc_ODE_prealloc(self,parameters,TimeVar,case,N,x0,t,error):
        errorcheck = 0
        
        
        while errorcheck == 0:
            FSP = lambda g,x: np.ravel(np.dot(At[0][g],np.array([x]).T))
            
            #ode15s.set_integrator('vode', method='bdf', order=15, nsteps=3000)
            At = self.Test_A(parameters,N,TimeVar,case,t)
            Jac = lambda x,t: At[t]
            ode15s = ode(FSP,jac=Jac)
            ode15s.set_integrator('dopri5',nsteps=3000)
            ode15s.set_initial_value(x0, t[0]) 
            solution = np.empty((len(t),3*N))
            solution[0,:] = x0[:,0]
        
            for g in range(1,len(t)):
                ode15s.set_initial_value(solution[g-1,:],t[g-1])
                                
                solution[g,:] = ode15s.integrate(g)            
            errorcheck = 1
            for i in range(0,len(t)):
                if 1-sum(solution[i]) > error:
                    errorcheck = 0
                    
            if errorcheck == 0:
                N = int(np.ceil(1.1*N))
                
                x0i = x0[0:3]
                x0 = np.zeros((1,3*N)).T
                x0[0:3] = x0i
                
            if errorcheck == 1:
                
                return solution,N
                
            if N > 1000:
                break


    #Not used
    
    def FSPCalc_ODE(self,parameters,TimeVar,case,N,x0,t,error):
        errorcheck = 0
        
        while errorcheck == 0:
        
            A = lambda t: self.Create_CME_TimeVar(parameters,N,TimeVar,case,t)
            
            FSP = lambda t,x: np.ravel(np.dot(A(t),np.array([x]).T))
            ode15s = ode(FSP)
            #ode15s.set_integrator('vode', method='bdf', order=15, nsteps=3000)
            ode15s.set_integrator('dop853')
            ode15s.set_initial_value(x0, t[0])
            solution = np.empty((len(t),3*N))
            solution[0,:] = x0[:,0]
            for i in range(1,len(t)):
                solution[i,:] = ode15s.integrate(t[i])[:,0]
            
            errorcheck = 1
            for i in range(0,len(t)):
                if 1-sum(solution[i]) > error:
                    errorcheck = 0
                    
            if errorcheck == 0:
                N = np.ceil(N*1.1)
                x0i = x0[0:3]
                x0 = np.zeros((1,3*N)).T
                x0[0:3] = x0i
                

            if errorcheck == 1:
                return solution
                
            if N > 1000:
                break

    def FSP_smooth_graph(self,solution,N): 
        
        RNAprob = np.zeros((1,N))    
        for i in range(0,len(RNAprob[0][:])-1):            
            RNAprob[0][i] = sum(solution[3*i:3*i+3]) 
            
        x = np.linspace(0,N-1,N)
        x_sm = np.array(x)
        y = RNAprob.flatten()
        y_sm = np.array(y)
        x_smooth = np.linspace(x_sm.min(), x_sm.max(), 200)
        y_smooth = spline(x, y, x_smooth)
        
        return x_smooth,y_smooth

    #Return probability 
    def FSP_Return_RNA_Prob(self,solution,N):
        RNAprob = np.zeros((1,N))    
        for i in range(0,len(RNAprob[0][:])-1):            
            RNAprob[0][i] = sum(solution[3*i:3*i+3]) 
            
        RNAprob = RNAprob.flatten()
        c = np.linspace(0,N-1,N)
        return c,RNAprob

            
    
    def TimeVar_Create(self,string):
        
        ##small parsing to prevent code execution inside the future execute call
        if ";" in string:
            return
        if "import" in string:
            return
        if "eval" in string:
            return
        if "exec" in string:
            return
        if "run" in string:
            return
        if "os." in string:
            return
        if ".py" in string:
            return
            
        string = string.lower()
        string = string.replace("^","**")
        string = string.replace("pi","np.pi")
        
        string = string.replace("sin","np.sin")
        string = string.replace("cos","np.cos")
        string = string.replace("tan","np.tan")
        string = string.replace("e^","np.exp")
        string = string.replace("exp","np.exp")
        string = string.replace("log","np.log")
        string = string.replace("ln","np.log")
        string = string.replace("sinh","np.sinh")
        string = string.replace("cosh","np.cosh")
        string = string.replace("tanh","np.tanh")        
        string = string.replace("arcsinh","np.arcsinh")
        string = string.replace("arccosh","np.arccosh")
        string = string.replace("arctanh","np.arctanh")
        string = string.replace("arcsin","np.arcsin")
        string = string.replace("arccos","np.arccos")
        string = string.replace("arctan","np.arctan")
        
        string = string.replace("log10","np.log10")
        string = string.replace("log2","np.log2")
        
        exec("TimeVar = lambda t:"+ string)
        
        
        return TimeVar
        
    #Ode model of 3 state gene with 2 transcription states and 1 off state
    def Gene_Model_A(self,t,k12,k23,k21,k32,kr2,kr3,gamma,beta):
        
        A = np.array([[-k12,k21,0,0], #in out of state 1
             [k12,-k21-k23,k32,0], #in out of state 2
             [0,k23,-k32,0],#in out of state 3
             [0,kr2,kr3,-gamma]] )
        return A
        
    def Gene_Model_A_TimeVar(self,t,k12,k23,k21,k32,kr2,kr3,gamma,beta,case,TimeVar):
        
        if case == 1:
            k12 = max(0,k12 +beta*TimeVar(t))
        if case == 2:
            k23 = max(0,k23 + beta*TimeVar(t))
        if case == 3:
            k21 = max(0,k21+ beta*TimeVar(t))
        if case == 4:
            k32 = max(0,k32 + beta*TimeVar(t))
        
        
        A = np.array([[-k12,k21,0,0], #in out of state 1
                 [k12,-k21-k23,k32,0], #in out of state 2
                 [0,k23,-k32,0],#in out of state 3
                 [0,kr2,kr3,-gamma]] )
                 

        return A
        
    def Gene_Model_Jac(self,x,t,k12,k23,k21,k32,kr2,kr3,gamma,beta,case,TimeVar):
        
        if case == 1:
            k12 = max(k12 +beta*TimeVar(t),0)
        if case == 2:
            k23 = max(0,k23 + beta*TimeVar(t))
        if case == 3:
            k21 = max(0,k21+ beta*TimeVar(t))
        if case == 4:
            k32 = max(0,k32 + beta*TimeVar(t))   
            
            
        A = np.array([[-k12,k21,0,0], #in out of state 1
             [k12,-k21-k23,k32,0], #in out of state 2
             [0,k23,-k32,0],#in out of state 3
             [0,kr2,kr3,-gamma]] )
      
        return A                
                
                
    def Gene_Model(self,x,t,k12,k23,k21,k32,kr2,kr3,gamma,beta):
            
        A = np.array([[-k12,k21,0,0], #in out of state 1
                 [k12,-k21-k23,k32,0], #in out of state 2
                 [0,k23,-k32,0],#in out of state 3
                 [0,kr2,kr3,-gamma]] )      
    #ode vector
        dxdt = np.dot(A,np.array([x]).T)
        dxdt = dxdt.flatten()
        return dxdt
        
    #Solver for the non time varying case
    def Gene_Dot(self,x0,t,parameters):

        
        solution = odeint(self.Gene_Model, x0,t,args = tuple(parameters))
        
        return solution
        
    #Solver for the time varying model
    def Gene_Dot_TimeVar(self,x0,t,parameters,case,TimeVar):
        
        Jac = self.Gene_Model_Jac
        solution = odeint(self.Gene_Model_TimeVar_compressed, x0,t,Dfun = Jac,args = tuple(parameters)+tuple([case])+tuple([TimeVar]))

        return solution
                
        
        
    
        
    def Gene_Model_TimeVar(self,x,t,k12,k23,k21,k32,kr2,kr3,gamma,beta,case,TimeVar):
        Tpass = TimeVar(t)
        if case ==1:
                k12 = max(k12 +beta*Tpass,0)
                S1,S2,S3,R = x #the 3 states and R (RNA)
    #ode vector
                dxdt = [k21*S2-k12*Tpass*S1,       #in out of state 1
                k12*TimeVar(t)*S1+k32*S3-k21*S2-k23*S2, #in out of state 2
                k23*S2-k32*S3,                          #in out of state 3
                kr2*S2+kr3*S3-gamma*R]                  #in out of RNA
                
                
        if case ==2:   
                k23 = max(0,k23 + beta*Tpass)
                
                S1,S2,S3,R = x #the 3 states and R (RNA)
    #ode vector
                dxdt = [k21*S2-k12*S1,                  #in out of state 1
                k12*S1+k32*S3-k21*S2-k23*Tpass*S2, #in out of state 2
                k23*TimeVar(t)*S2-k32*S3,               #in out of state 3
                kr2*S2+kr3*S3-gamma*R]                  #in out of RNA  
                
        if case ==3:
                k21 = max(0,k21+ beta*Tpass)
                S1,S2,S3,R = x 
            #ode vector
                dxdt = [k21*Tpass*S2-k12*S1,       #in out of state 1
                k12*S1+k32*S3-k21*Tpass*S2-k23*S2, #in out of state 2
                k23*S2-k32*S3,                          #in out of state 3
                kr2*S2+kr3*S3-gamma*R]                  #in out of RNA  
                
        if case ==4:
                k32 = max(0,k32 + beta*Tpass)
                S1,S2,S3,R = x 
                dxdt = [k21*S2-k12*S1,                  #in out of state 1
                k12*S1+k32*Tpass*S3-k21*S2-k23*S2, #in out of state 2
                k23*S2-k32*Tpass*S3,               #in out of state 3
                kr2*S2+kr3*S3-gamma*R]                  #in out of RNA  
        
        return dxdt
    




        
    def Gene_Model_TimeVar_compressed(self,x,t,k12,k23,k21,k32,kr2,kr3,gamma,beta,case,TimeVar):
        Tpass = TimeVar(t)
        if case ==1:
                k12 = max(k12 +beta*Tpass,0)
                
                
        if case ==2:   
                k23 = max(0,k23 + beta*Tpass)
                

        if case ==3:
                k21 = max(0,k21+ beta*Tpass)
                

        if case ==4:
                k32 = max(0,k32 + beta*Tpass)
     #in out of RNA  
                
        
        A = np.array([[-k12,k21,0,0], #in out of state 1
         [k12,-k21-k23,k32,0], #in out of state 2
         [0,k23,-k32,0],#in out of state 3
         [0,kr2,kr3,-gamma]] )      
#ode vector
        dxdt = np.dot(A,np.array([x]).T)
        dxdt = dxdt.flatten()

        
        return dxdt






    
    #Create a function to import our csv file of the mRNA histogram
    def Import_mRNA(self,mRNA_Hist_file):
     #Use loadtext from numpy to put it directly in an array
        mRNA_arr = np.loadtxt(mRNA_Hist_file, delimiter = ',') 
        return mRNA_arr
        
    def mRNA_Hist(self,mRNA_arr):
        x=1
        
    def Solve_Moments(self,x0,t,parameters):

        
        solution = odeint(self.Second_Moments, x0,t,args = tuple(parameters))
        
        return solution
        
    #Solver for the time varying model
    def Solve_Moments_TimeVar(self,x0,t,parameters,case,TimeVar):
        

        solution = odeint(self.Second_Moments_TimeVar, x0,t,args = tuple(parameters)+tuple([case])+tuple([TimeVar]))

        return solution
                
    def Second_Moments(self,x,t,k12,k23,k21,k32,kr2,kr3,gamma,beta):
        
        
                
        secondmoments = np.array([[-k12,        k21,    0,      0,      0,                0,          0,            0,              0,                0,                  0,      0,            0,        0],
                        [ k12, -k21 - k23,  k32,      0,      0,                0,          0,            0,              0,                0,                  0,      0,            0,        0],
                        [   0,        k23, -k32,      0,      0,                0,          0,            0,              0,                0,                  0,      0,            0,        0],
                        [   0,        kr2,  kr3, -gamma,      0,                0,          0,            0,              0,                0,                  0,      0,            0,        0],
                        [ k12,        k21,    0,      0, -2*k12,              k21,          0,            0,              0,                0,                  0,      0,            0,        0],
                        [-k12,       -k21,    0,      0,    k12, -k12 - k21 - k23,        k32,            0,            k21,                0,                  0,      0,            0,        0],
                        [   0,          0,    0,      0,      0,              k23, -k12 - k32,            0,              0,              k21,                  0,      0,            0,        0],
                        [   0,          0,    0,      0,      0,              kr2,        kr3, -gamma - k12,              0,                0,                k21,      0,            0,        0],
                        [ k12,  k21 + k23,  k32,      0,      0,              k12,          0,            0, -2*k21 - 2*k23,              k32,                  0,      0,            0,        0],
                        [   0,       -k23, -k32,      0,      0,                0,        k12,            0,            k23, -k21 - k23 - k32,                  0,    k32,            0,        0],
                        [   0,          0,    0,      0,      0,                0,          0,          k12,            kr2,              kr3, -gamma - k21 - k23,      0,          k32,        0],
                        [   0,        k23,  k32,      0,      0,                0,          0,            0,              0,              k23,                  0, -2*k32,            0,        0],
                        [   0,          0,    0,      0,      0,                0,          0,            0,              0,                0,                k23,    kr3, -gamma - k32,        0],
                        [   0,        kr2,  kr3,  gamma,      0,                0,          0,            0,              0,                0,                kr2,      0,          kr3, -2*gamma]])
                        
        dxdt = np.dot(secondmoments,np.array([x]).T)
        dxdt = dxdt.flatten()
        return dxdt    
    
    def Second_Moments_TimeVar(self,x,t,k12,k23,k21,k32,kr2,kr3,gamma,beta,case,TimeVar):
        
        
        if case !=0:
                Tpass = TimeVar(t)
                
        if case ==1:
                k12 = max(k12 +beta*Tpass,0)
                
                
        if case ==2:   
                k23 = max(0,k23 + beta*Tpass)
                

        if case ==3:
                k21 = max(0,k21+ beta*Tpass)
                

        if case ==4:
                k32 = max(0,k32 + beta*Tpass)
                
        secondmoments = np.array([[-k12,        k21,    0,      0,      0,                0,          0,            0,              0,                0,                  0,      0,            0,        0],
                        [ k12, -k21 - k23,  k32,      0,      0,                0,          0,            0,              0,                0,                  0,      0,            0,        0],
                        [   0,        k23, -k32,      0,      0,                0,          0,            0,              0,                0,                  0,      0,            0,        0],
                        [   0,        kr2,  kr3, -gamma,      0,                0,          0,            0,              0,                0,                  0,      0,            0,        0],
                        [ k12,        k21,    0,      0, -2*k12,              k21,          0,            0,              0,                0,                  0,      0,            0,        0],
                        [-k12,       -k21,    0,      0,    k12, -k12 - k21 - k23,        k32,            0,            k21,                0,                  0,      0,            0,        0],
                        [   0,          0,    0,      0,      0,              k23, -k12 - k32,            0,              0,              k21,                  0,      0,            0,        0],
                        [   0,          0,    0,      0,      0,              kr2,        kr3, -gamma - k12,              0,                0,                k21,      0,            0,        0],
                        [ k12,  k21 + k23,  k32,      0,      0,              k12,          0,            0, -2*k21 - 2*k23,              k32,                  0,      0,            0,        0],
                        [   0,       -k23, -k32,      0,      0,                0,        k12,            0,            k23, -k21 - k23 - k32,                  0,    k32,            0,        0],
                        [   0,          0,    0,      0,      0,                0,          0,          k12,            kr2,              kr3, -gamma - k21 - k23,      0,          k32,        0],
                        [   0,        k23,  k32,      0,      0,                0,          0,            0,              0,              k23,                  0, -2*k32,            0,        0],
                        [   0,          0,    0,      0,      0,                0,          0,            0,              0,                0,                k23,    kr3, -gamma - k32,        0],
                        [   0,        kr2,  kr3,  gamma,      0,                0,          0,            0,              0,                0,                kr2,      0,          kr3, -2*gamma]])
                        
        dxdt = np.dot(secondmoments,np.array([x]).T)
        dxdt = dxdt.flatten()
        return dxdt
        
    
    
    
    def StoichFun(self):  #Set up a function to return the Stoichiometry matrix
        
        #             1>2 2>3 2>1 3>2 2>R 3>R R>degrade
        S = np.array([[-1, 0, 1,  0,  0,   0,  0],  #S1
                     [  1,-1,-1,  1,  0,   0,  0],  #S2
                     [  0, 1, 0, -1,  0,   0,  0],  #S3
                     [  0, 0, 0,  0,  1,   1, -1]]) #R
        
        return S
    
    
    
    
    
    
    #Function to set up Propensity based on defined parameters
    def PropensityFun(self,t,x,parameters,case,TimeVar):
    
        k12 = parameters[0]
        k23 = parameters[1]
        k21 = parameters[2]
        k32 = parameters[3]
        kr2 = parameters[4]
        kr3 = parameters[5]
        beta = parameters[7]
        gamma = parameters[6]
        Tpass = t.item(0,0)
        
        if case == 1:
            k12 = max(parameters[0] +beta*TimeVar(Tpass),0)
        if case == 2:
            k23 = max(0,parameters[1] + beta*TimeVar(Tpass))
        if case == 3:
            k21 = max(0,parameters[2]+ beta*TimeVar(Tpass))
        if case == 4:
            k32 = max(0,parameters[3] + beta*TimeVar(Tpass))
        
        #Propensitys 7 reactions 1>2 2>3 2>1 3>2 Kr2 KR3 R>degrade
        prop = np.array([
                    [k12*x[0]],
                    [k23*x[1]],
                    [k21*x[1]],
                    [k32*x[2]],
                    [kr2*x[1]],
                    [kr3*x[2]],
                    [gamma*x[3]]])   
        return prop
    
        
    def SSA(self,x0,Output_Times,parameters,case,TimeVar):
        t = np.array([[float(0)]]) #Inital Time
        x = np.copy(x0)  #Copy the x0 so as not to change x0 over the SSA
        ip = 0 #Index Pointer for Output_Times
        S = self.StoichFun()#Get the S matrix
        nt = (4,len(Output_Times)) 
        X_Array = np.zeros(nt) #Preallocate X_Array
        
    
        while t < Output_Times[-1:]: #While t < final time
            
            w = self.PropensityFun(t,x,parameters,case,TimeVar) #Update propensities
    
            w0 = sum(w) #get the sum
            
           # tau = -np.log(np.random.random(1))/w0 #randomly generate next rxn time
            
            tau = np.random.exponential(1.0/w0)
            
            r2 = w0*np.random.random(1) #randomly generate which reaction happens
            u = 1
            
            while sum(w[0:u,0]) < r2:    #Check to see which reaction happened      
                u = u + 1
            
                
            t = t+ tau  #update time
            while ip < len(Output_Times) and t > Output_Times[ip]: #Record the curent X if t > current output times index
                X_Array[:,ip] = x[:,0]
                ip = ip+1
            
            
            if ip == len(Output_Times):
                return X_Array
                
                
            x[:,0] = x[:,0]+S[:,u-1]
            #Update the X matrix based on uth reaction
        
    #Take the average of Nruns of SSA
    def SSA_u(self,x0,t,parameters,Nruns,case,TimeVar):
        nt = (len(t),Nruns) #Preallocate the matrices
        S1 = np.zeros(nt)
        S2 = np.zeros(nt)
        S3 = np.zeros(nt)
        R = np.zeros(nt)
           
        for i in range(0,Nruns):
            #For the range of Nruns run an SSA and record it
            X_Array = self.SSA(x0,t,parameters,case,TimeVar)
            S1[:,i] = X_Array[0,:]
            S2[:,i] = X_Array[1,:]
            S3[:,i] = X_Array[2,:]
            R[:,i] = X_Array[3,:]
        #take the average across the Nrun dimension of the recorded runs
        S1_u = np.average(S1,1)
        S2_u = np.average(S2,1)
        S3_u = np.average(S3,1)
        R_u = np.average(R,1)
        averages = np.array([S1_u,S2_u,S3_u,R_u])
        #Plot the figures
        
        return averages    
        
        
        
    def SSA_u_compressed(self,x0,t,parameters,Nruns,case,TimeVar):
        solution_mat = np.empty((len(x0),len(t),Nruns))
        for i in range(Nruns):
            solution_mat[:,:,i] = self.SSA(x0,t,parameters,case,TimeVar)
        
        
    #Create a function to give the difference of the solution data and the experimental data
    def DifferenceODE(self,Parameters,case,Output_Times,x0,Data_Set,TimeVar,standevData_set,meanData_set):
       
        if case in [1,2,3,4]: #if its a time varying case use the timevarying function
            
            solution = self.Gene_Dot_TimeVar(x0,Output_Times,Parameters,case,TimeVar) 
        else: #otherwise dont use the time varying case
            solution = self.Gene_Dot(x0,Output_Times,Parameters)
            #reshape our for easier computing
        
        mRNA_solution = np.reshape(solution[:,3],(1,len(Data_Set[0,:])))
        #mRNA_solution = solution[:,3].T
        sseornot = sum(standevData_set[0,:])
        print len(Data_Set[0,:])
        
        #intialize the parameters 
        diff = 0
        #for the length of the data set caluclate the average/std over time for all cells
        for i in range(0,len(Data_Set[0,:])):
            # if its not the first row (all zeros) calculate the standard sum of errors squared
            if i > 0:
                if sseornot==0:
                    diff = diff + (meanData_set[0,i] - mRNA_solution[0,i])**2 
                else:
                    diff = diff + ((meanData_set[0,i] - mRNA_solution[0,i])/(standevData_set[0,i]/len(Data_Set[0,:])))**2 
                    #Change to standard error about the mean
        
        
        
        return diff


    def Difference_Second_Moments(self,Parameters,case,Output_Times,x0,Data_Set,TimeVar,standevData_set,meanData_set):
    
        if case in [1,2,3,4]: #if its a time varying case use the timevarying function
            
            solution = self.Solve_Moments_TimeVar(x0,Output_Times,Parameters,case,TimeVar)
        else: #otherwise dont use the time varying case
            solution = self.Solve_Moments(x0,Output_Times,Parameters)
            #reshape our for easier computing
        
        #Variance = solution[:,13]        
        
        mRNA_mean = np.reshape(solution[:,3],(1,len(Data_Set[0,:])))
        mRNA_var = np.reshape(solution[:,13],(1,len(Data_Set[0,:])))
        #mRNA_solution = solution[:,3].T
        sseornot = sum(standevData_set[0,:])

        
        #intialize the parameters 
        diff = 0
        
        #for the length of the data set caluclate the average/std over time for all cells
        for i in range(0,len(Data_Set[0,:])):
            # if its not the first row (all zeros) calculate the standard sum of errors squared
            if i > 0:
               diff = diff + -( -(len(Data_Set[:,0])-1)/2.0*np.log(mRNA_var[0,i]) - (len(Data_Set[:,0])-1)*standevData_set[0,i]**2/(2*mRNA_var[0,i]))
                    
        
        
        
        return diff
    
    def MinDifferenceODE(self,case,Output_Times,x0,Data_Set,Parameter_Guess,TimeVar,standevData_set,meanData_set,Bounds,fminmethod):
        #Redefine our difference function as an autonomous function of paramters (P)
        
        objectivefun = lambda p: self.DifferenceODE(p,case,Output_Times,x0,Data_Set,TimeVar,standevData_set,meanData_set)
        bounds = self.Get_Par_Bounds(case)
        bounds = np.float64(bounds)
        best_parameters = minimize(objectivefun,Parameter_Guess,method=fminmethod,bounds = Bounds,options ={'maxiter':500}) #use the built in optimiztion of sci py
        return best_parameters #return the best parameters
        

    def MinDifferenceFSP(self,case,Bounds,fminmethod,Function,Parameter_Guess):
        #Redefine our difference function as an autonomous function of paramters (P)
        
        objectivefun = Function
        bounds = self.Get_Par_Bounds(case)
        bounds = np.float64(bounds)
        
        best_parameters = minimize(objectivefun,Parameter_Guess,method=fminmethod,bounds = Bounds,options ={'maxiter':500}) #use the built in optimiztion of sci py
        return best_parameters #return the best parameters
        

        
    def Second_Moment_fit(self,case,Output_Times,x0,Data_Set,Parameter_Guess,TimeVar,standevData_set,meanData_set,Bound,fminmethod):
        
        objectivefun = lambda p: self.DifferenceODE(p,case,Output_Times,x0,Data_Set,TimeVar,standevData_set,meanData_set)
        bounds = self.Get_Par_Bounds(case)
        bounds = np.float64(bounds)
        best_parameters = minimize(objectivefun,Parameter_Guess,method=fminmethod,bounds = Bounds,options ={'maxiter':500}) #use the built in optimiztion of sci py
        return best_parameters #return the best parameters
                

    def make_binary(self,x,xmin,xmax,N):
        #convert a number to binary
        y = float(x-xmin)/(xmax-xmin)
        binary = np.zeros(N)
        for i in range(0,N):
            if y >= .5:
                binary[i] = 1
                y = y - .5
            y = y*2
        return binary
    
    def make_binaryvec(self,x,xmin,xmax,N,M):
        #encode a binary matrix for faster met-haste searches / SA searches
        binaryvec = np.zeros(M*N)
        for i in range(0,M):
            binaryvec[(i)*N:(i+1)*N] = self.make_binary(x[i],xmin[i],xmax[i],N)
        return binaryvec
        
    def make_floatvec(self,binary,xmin,xmax,N,M):
        #Function to converty binary vector back to float
        fvec = np.zeros(M)
        bincon = np.zeros(N)
    
        for i in range(0,N):
            bincon[i] = N-(i+1)
            
        for i in range(0,M):
                floatvec = binary[i*N:(i+1)*N]
                fvec[i] = xmin[i] + sum((2**bincon*floatvec))/(2**N-1)*(xmax[i]-xmin[i])
        return fvec
        
        
            
    def Get_Par_Bounds(self,case):
        ## Function to return bounds for parameter minimizations by model case
        from numpy import zeros as zr
        Bounds =  zr((2,8))
        Bounds[0,0:7] = .1
        Bounds[1,0:3] = 30
        Bounds[1,3:7] = 30
        Bounds[0,7] =-10
        Bounds[1,7] = 10
        if case == 1:
            Bounds[0,0] = -2
            Bounds[1,0] = 10
            Bounds[0,7] = .1
            Bounds[1,7] = 10
        if case == 3:
            Bounds[0,1] = .1
            Bounds[1,1] = 10
            Bounds[0,7] = -2
            Bounds[1,7] = .1
        if case == 2:
            Bounds[0,2] = -2
            Bounds[1,2] = 10
            Bounds[0,7] = -2
            Bounds[1,7] = .1
        if case == 4:
            Bounds[0,3] = .1
            Bounds[1,3] = 10
            Bounds[0,7] = -2
            Bounds[1,7] = .1
        return Bounds
          
    
    
    def MetroplisHastings(self,Function,Init_Parameters,Bounds,N_Burn,N_Chain,N_Thin,Mut_rate):
        x0 = np.copy(Init_Parameters)
        x = np.copy(x0)#self.make_binaryvec(x0,Bounds[0,:],Bounds[1,:],N,M) #convert the intial parameter guess to binary
        xbest = np.copy(x) #intialize the xbest storage variable
        f = Function #define the function
        #xfloat = self.make_floatvec(x,Bounds[0,:],Bounds[1,:],N,M) #intialize xfloat
        fx = f(x) #intialize fx
        fbest = np.copy(fx) #initalize fbest
        funchain = np.array([]) #intialize output variables
        parchain = Init_Parameters
    
        for i in range(-N_Burn,N_Chain): #for the length of Nchain, record
            
            for j in range(0,N_Thin):
    
                xnew = np.copy(x) #copy the binary vector
               
                #mutate random bases 0 > 1 or 1 > 0
#                mut_array = np.where(np.random.random((N*M)) <Mut_rate)
#                for k in range(0,len(mut_array[0])):
#                    xnew[mut_array[0][k]] = np.abs(xnew[mut_array[0][k]]-1)
#                
#               
#                xfloat = self.make_floatvec(xnew,Bounds[0,:],Bounds[1,:],N,M)
                for k in range(0,len(xnew)):
                    if np.random.random() < Mut_rate:
                        xnew[k] = (Bounds[1,k] - Bounds[0,k])*np.random.random()+Bounds[0,k]

                       
                        
                fnew = f(xnew)
                
                if fnew < fx or np.random.random(1) <= np.exp(fx - fnew):
                    x = np.copy(xnew)
                    fx = np.copy(fnew)
                    
                    if fnew < fbest:
                        fbest = np.copy(fnew)
                        xbest = np.copy(xnew)
                      
                if i >= 0:
                    
                    if j == N_Thin-1:
                    
                        funchain = np.append(funchain,fx)
                        parchain = np.vstack([parchain, x])
                        
                        
        return funchain, parchain, fbest, xbest
        
        
    def SimulatedAnnealing(self,Function,Init_Parameters,Bounds,mut_rate):
        x0 = np.copy(Init_Parameters)
        M = len(x0)
        N = 64
        funchain = np.array([])
        parchain = Init_Parameters    
        x = self.make_binaryvec(x0,Bounds[0,:],Bounds[1,:],N,M)
        xfloat = np.copy(x)
        xbest = np.copy(x)
        f = Function
        
        xfloat = self.make_floatvec(x,Bounds[0,:],Bounds[1,:],N,M)
        
        fx = f(xfloat)
        fbest = np.copy(fx)   
        T=100
        Tmin = 1e-10
        dT = .001
        while T>Tmin:
            xnew = np.copy(x)
            mut_array = np.where(np.random.random((N*M)) <mut_rate)
            for k in range(0,len(mut_array[0])):
                if xnew[mut_array[0][k]] == 1:
                    xnew[mut_array[0][k]] = 0
                else:
                    
                    j = mut_array[0][k]
                    
                    while xnew[j] == 0:
                        xnew[j] = 1
                        j += 1
                        if xnew[j] == 1 and j < 128:
                            xnew[j]= 0
                
                        
                
            xfloat = self.make_floatvec(xnew,Bounds[0,:],Bounds[1,:],N,M)
            fnew = f(xfloat)
            
            if fnew <=fx or np.log(np.random.random(1)) <= ((fx-fnew)/T):
               
                x = np.copy(xnew)
                fx = np.copy(fnew)
    
                
                if fnew < fbest:
                    fbest = np.copy(fnew)
    
                    xbest = np.copy(xnew)
                           
                funchain = np.append(funchain,fx)
                parchain = np.vstack([parchain, xfloat])        
            T = T*(1-dT)
            
        floatbest = self.make_floatvec(xbest,Bounds[0,:],Bounds[1,:],N,M)
        
        return floatbest,fbest,funchain,parchain

class GUI(Frame): #define the GUI class


    ##########################################################################
    #initialize class variables and other UI
    ##########################################################################
    def mouse_location(self,event):
        #continiously update the mouse coordinates from motion across canvas
        self.mouse_x, self.mouse_y = event.x,event.y
        
        
        
    def __init__(self, parent): #define the intial setup of the GUI
        
        Frame.__init__(self, parent)   
        self.GeneModel = Gene_Model_Class.GeneModel()
        self.parent = parent
        
        #custom font and colormap
        self.customFont = tkFont.Font(size=10)
        self.colormap = ['blue','red','green','magenta','cyan', 'yellow', 'black','#00FF00']
        self.customcolormap = ['blue','red','green','magenta','cyan', 'yellow', 'black','#00FF00']
        self.Open_Fig() #open up a figure
        self.bounds = self.GeneModel.Get_Par_Bounds(2) #initalize bounds
        self.MH_options_open = False
        #self.progress = Progressbar(LoadingandFigFrame,orient=HORIZONTAL, length=241, mode='determinate')
       # self.progress.grid(row=0,column=0,padx=2,pady=2)
    
        self.holdon = IntVar() #hold on variable for plotting multiple
        self.holdon.set(0)
        self.x = 0
        self.y = 0
        self.mouse_x = 1
        self.mouse_y = 1
        self.count = 0
        self.parent.bind('<Enter>',self.start_poll)
        self.parent.bind('<Leave>',self.stop_poll)
        self.offscreen = True
        self.menuopen = False
        self.colormap_select_var = IntVar() #colormap select variable 0 - 4
        self.colormap_select_var.set(0)
        
        self.parent.bind('<Motion>',self.mouse_location)
        self.initUI() 
        
    
###############################################################################
# MATH METHODS
# math methods defined at the top of the GUI
###############################################################################


    def SSA_multi(self,runs):
        Output_Times =self.Output_Times_ssa_pass
        x0 = self.x0_ssa_pass
        case = self.case_ssa_pass
        parameters = self.parameters_ssa_pass
        t = np.array([[float(0)]]) #Inital Time
        x = np.copy(x0)  #Copy the x0 so as not to change x0 over the SSA
        TimeVar= self.TimeVar_ssa_pass
        ip = 0 #Index Pointer for Output_Times
        S = self.GeneModel.StoichFun()#Get the S matrix
        nt = (4,len(Output_Times)) 
        X_Array = np.zeros(nt) #Preallocate X_Array
        
    
        while t < Output_Times[-1:]: #While t < final time
            
            w = self.GeneModel.PropensityFun(t,x,parameters,case,TimeVar) #Update propensities
    
            w0 = sum(w) #get the sum
            
           # tau = -np.log(np.random.random(1))/w0 #randomly generate next rxn time
            
            tau = np.random.exponential(1.0/w0)
     
           
            
            r2 = w0*np.random.random(1) #randomly generate which reaction happens
            u = 1
            
            while sum(w[0:u,0]) < r2:    #Check to see which reaction happened      
                u = u + 1
            
                
            t = t+ tau  #update time
            while ip < len(Output_Times) and t > Output_Times[ip]: #Record the curent X if t > current output times index
    
                X_Array[:,ip] = x[:,0]
                ip = ip+1
            
            
            if ip == len(Output_Times):
                return X_Array
                
            
            x[:,0] = x[:,0]+S[:,u-1]


    ########################################################################
    # Average SSA runs
    # run the average SSA runs in the top level of the GUI for speed
    # and progress bar updating
    ########################################################################
    def SSA_u(self,x0,t,parameters,Nruns,case,TimeVar):
        
        #define a popup progress bar to take focus so you dont have to update all the GUI
        if Nruns > 200:
            self.methasteprog_window = Toplevel()
            self.methasteprog_window.title('SSA Progress')
            self.methasteprog_window.iconbitmap('qbio_icon_9nD_icon.ico')
            self.methasteprog_window.focus_force()
            self.methasteprog_window.grab_set()
            methasteprog_frame = Frame(self.methasteprog_window)
            self.methasteprog = Progressbar(methasteprog_frame,orient=HORIZONTAL, length=400, mode='determinate')
            methasteprog_frame.pack()
            
            self.methasteprog_window.resizable(width=False,height=False)
            starttime = time.time()
            now = starttime
            evalpersec = 1
            evalmax = 1
            predvar = 0
            predicttime = []
            self.methasteprog.grid(row=0,column=0,columnspan=4)
            self.prog_clock = Label(methasteprog_frame,text='')
            self.prog_clock.grid(row=1,column=1,pady=2,padx=2,sticky=W)
            prog_label = Label(methasteprog_frame,text='Run Time:')
            pred_label = Label(methasteprog_frame,text='Prediction:')
            prog_label.grid(row=1,column=0,pady=2,padx=2,sticky=E)
            pred_label.grid(row=1,column=2,pady=2,padx=2,sticky=E)
            
            self.pred_clock = Label(methasteprog_frame,text='')
            self.pred_clock.grid(row=1,column=3,pady=2,padx=2,sticky=W)
            
            
            self.methasteprog.grid(row=0,column=0)
            self.methasteprog["maximum"]=Nruns #set max of prog bar
            self.methasteprog.grab_set()
        
        nt = (len(t),Nruns) #Preallocate the matrices
        S1 = np.zeros(nt)
        S2 = np.zeros(nt)
        S3 = np.zeros(nt)
        R = np.zeros(nt)
        
        self.x0_ssa_pass =x0
        self.parameters_ssa_pass = parameters
        self.Output_Times_ssa_pass = t
        self.case_ssa_pass = case
        self.TimeVar_ssa_pass = TimeVar
        '''
        pool = ThreadPool(10) 
        results = pool.map(self.SSA_multi,range(Nruns))
        
        for i in range(0,len(results)):
            S1[:,i] = results[i][:][0]
            S2[:,i] =results[i][:][1]
            S3[:,i] = results[i][:][2]
            R[:,i] = results[i][:][3]
        '''
        for i in range(0,Nruns):
            #For the range of Nruns run an SSA and record it
            X_Array = self.GeneModel.SSA(x0,t,parameters,case,TimeVar)
            S1[:,i] = X_Array[0,:]
            S2[:,i] = X_Array[1,:]
            S3[:,i] = X_Array[2,:]
            R[:,i] = X_Array[3,:]
            if Nruns > 200:
                self.methasteprog.step(amount=1)
                nowold = now
                now = int(np.floor(time.time()-starttime))
                
                
                if nowold == now:
                    evalpersec= evalpersec + 1     
                    
                else:
                    if len(predicttime) < 5:
                        if len(predicttime) == 1:
                            predicttime = []
                        predicttime.append(int((Nruns)/evalpersec))
                    
                        predvar= sum(predicttime)/len(predicttime)
                    else:
                        predicttime.append(int((Nruns)/evalpersec))
                        predicttime = predicttime[1:]
                        predvar= sum(predicttime)/len(predicttime)
                    evalpersec =1
                
                m, s = divmod(now, 60)
                h, m = divmod(m, 60)
                curr = "%d:%02d:%02d" % (h, m, s)
                self.prog_clock.configure(text=curr)
                
                m, s = divmod(predvar, 60)
                h, m = divmod(m, 60)
                pred = "%d:%02d:%02d" % (h, m, s)            
                self.pred_clock.configure(text=pred)
                self.methasteprog_window.update()   #update progressbar
            
        #take the average across the Nrun dimension of the recorded runs
        S1_u = np.average(S1,1)
        S2_u = np.average(S2,1)
        S3_u = np.average(S3,1)
        R_u = np.average(R,1)
        averages = np.array([S1_u,S2_u,S3_u,R_u])
        
        if Nruns > 200:
            self.methasteprog_window.destroy()    #destroy the progress bar window
            self.focus_force() #return focus to the main GUI
            
        return averages            

    ########################################################################
    # Gene Model integration, case 0 and TimeVariant
    # uses odeint (scipy) default
    ########################################################################

    def Gene_Dot(self,x0,t,parameters):
        solution = odeint(self.Gene_Model, x0,t,args = tuple(parameters))
        return solution
        
    def Gene_Dot_TimeVar(self,x0,t,parameters,case,TimeVar):
        solution = odeint(self.Gene_Model_TimeVar, x0,t, args = tuple(parameters)+tuple([case])+tuple([TimeVar]))
        return solution
 
    ########################################################################
    # Difference function to minimize
    # defined as calling the integration of the model and then subtracting
    # SSE or SSE/std if its a multi line (Hist) to fit
    ########################################################################
 
    def DifferenceODE(self,Parameters,case,Output_Times,x0,Data_Set,TimeVar,standevData_set,meanData_set):
        if case in [1,2,3,4]: #if its a time varying case use the timevarying function          
            solution = self.Gene_Dot_TimeVar(x0,Output_Times,Parameters,case,TimeVar) 
        else: #otherwise dont use the time varying case
            solution = self.Gene_Dot(x0,Output_Times,Parameters)
            #reshape our for easier computing      
        mRNA_solution = np.reshape(solution[:,3],(1,len(Data_Set[0,:])))
        
        sseornot = sum(standevData_set[0,:])     #check the file to fit
        #if the file is not a matrix histogram there is no standard dev, so dont use it
        
        diff = 0
        #for the length of the data set caluclate the average/std over time for all cells
        for i in range(0,len(Data_Set[0,:])):
            # if its not the first row (all zeros) calculate the standard sum of errors squared
            if i > 0:
                if sseornot==0:
                    diff = diff + (meanData_set[0,i] - mRNA_solution[0,i])**2 
                else:
                    diff = diff + ((meanData_set[0,i] - mRNA_solution[0,i])/standevData_set[0,i])**2           
        return diff
    
    ########################################################################
    # Define the Gene Models both case 0 and Time Varying
    # case 0 uses A*x notation
    # Time varying uses a more laid out notation
    ########################################################################
    
    def Gene_Model(self,x,t,k12,k23,k21,k32,kr2,kr3,gamma,beta):
            
        A = np.array([[-k12,k21,0,0], #in out of state 1
                 [k12,-k21-k23,k32,0], #in out of state 2
                 [0,k23,-k32,0],#in out of state 3
                 [0,kr2,kr3,-gamma]] )      
    #ode vector
        dxdt = np.dot(A,np.array([x]).T)
        dxdt = dxdt.flatten()
        return dxdt
    
    ########################################################################
    # Define the Gene Models both case 0 and Time Varying
    # case 0 uses A*x notation
    # Time varying uses a more laid out notation
    ########################################################################
         
    def Gene_Model_TimeVar(self,x,t,k12,k23,k21,k32,kr2,kr3,gamma,beta,case,TimeVar):
        Tpass = TimeVar(t) #save computing only call lambda function once
    
        if case ==1:
                k12 = max(k12 +beta*Tpass,0)
                S1,S2,S3,R = x #the 3 states and R (RNA)
    #ode vector
                dxdt = [k21*S2-k12*Tpass*S1,       #in out of state 1
                k12*TimeVar(t)*S1+k32*S3-k21*S2-k23*S2, #in out of state 2
                k23*S2-k32*S3,                          #in out of state 3
                kr2*S2+kr3*S3-gamma*R]                  #in out of RNA
                
                
        if case ==2:   
                k23 = max(0,k23 + beta*Tpass)
                
                S1,S2,S3,R = x #the 3 states and R (RNA)
    #ode vector
                dxdt = [k21*S2-k12*S1,                  #in out of state 1
                k12*S1+k32*S3-k21*S2-k23*Tpass*S2, #in out of state 2
                k23*TimeVar(t)*S2-k32*S3,               #in out of state 3
                kr2*S2+kr3*S3-gamma*R]                  #in out of RNA  
                
        if case ==3:
                k21 = max(0,k21+ beta*Tpass)
                S1,S2,S3,R = x 
            #ode vector
                dxdt = [k21*Tpass*S2-k12*S1,       #in out of state 1
                k12*S1+k32*S3-k21*Tpass*S2-k23*S2, #in out of state 2
                k23*S2-k32*S3,                          #in out of state 3
                kr2*S2+kr3*S3-gamma*R]                  #in out of RNA  
                
        if case ==4:
                k32 = max(0,k32 + beta*Tpass)
                S1,S2,S3,R = x 
                dxdt = [k21*S2-k12*S1,                  #in out of state 1
                k12*S1+k32*Tpass*S3-k21*S2-k23*S2, #in out of state 2
                k23*S2-k32*Tpass*S3,               #in out of state 3
                kr2*S2+kr3*S3-gamma*R]                  #in out of RNA  
        
        
        return dxdt
    
    
    ########################################################################
    # Define the MetHaste Search 
    # Defined in the top level of the GUI for speed
    # potential threading to be added
    ########################################################################            
            
    def MetroplisHastings(self,Function,Init_Parameters,Bounds,N_Burn,N_Chain,N_Thin,Mut_rate):
        #create a focus set progress bar to only update that instead of all the GUI
        self.methasteprog_window = Toplevel()
        self.methasteprog_window.iconbitmap('qbio_icon_9nD_icon.ico')
        self.methasteprog_window.title('Met-Haste Progress')
        self.methasteprog_window.focus_force()
        
        
        methasteprog_frame = Frame(self.methasteprog_window)
        self.methasteprog = Progressbar(methasteprog_frame,orient=HORIZONTAL, length=320, mode='determinate')
        methasteprog_frame.pack()
        self.methasteprog_window.resizable(width=False,height=False)
        starttime = time.time()
        now = starttime
        evalpersec = 1
        evalmax = 1
        predvar = 0
        predicttime = []
        self.methasteprog.grid(row=0,column=0,columnspan=4)
        self.prog_clock = Label(methasteprog_frame,text='')
        self.prog_clock.grid(row=1,column=1,pady=2,padx=2,sticky=W)
        prog_label = Label(methasteprog_frame,text='Run Time:')
        pred_label = Label(methasteprog_frame,text='Prediction:')
        prog_label.grid(row=1,column=0,pady=2,padx=2,sticky=E)
        pred_label.grid(row=1,column=2,pady=2,padx=2,sticky=E)
        
        self.pred_clock = Label(methasteprog_frame,text='')
        self.pred_clock.grid(row=1,column=3,pady=2,padx=2,sticky=W)
        
        self.methasteprog["maximum"]=N_Chain+N_Burn
        self.methasteprog.grab_set()
        
        x0 = np.copy(Init_Parameters)

        x = np.copy(x0)#self.make_binaryvec(x0,Bounds[0,:],Bounds[1,:],N,M) #convert the intial parameter guess to binary
        xbest = np.copy(x) #intialize the xbest storage variable
        f = Function #define the function
        #xfloat = self.make_floatvec(x,Bounds[0,:],Bounds[1,:],N,M) #intialize xfloat
        fx = f(x) #intialize fx
        fbest = np.copy(fx) #initalize fbest
        funchain = np.array([]) #intialize output variables
        parchain = Init_Parameters    
        waver_param = np.ones((1,8))
        accept_param = np.zeros((1,8))
        tot_param = np.zeros((1,8))
        waver_max = np.ones((1,8))
        for i in range(0,len(waver_max[0])):
            waver_max[0][i] = np.sqrt(abs(Bounds[1,i]-Bounds[0,i]))
        waver_min = np.ones((1,8))
        for i in range(0,len(waver_max[0])):
            waver_min[0][i] = abs(Bounds[1,i]-Bounds[0,i])*.0005
            

        for self.methastecounter in range(-N_Burn,N_Chain): #for the length of Nchain, record            
            for j in range(0,N_Thin):   #for len of N_thin
                xnew = np.copy(x) 
                
                '''
                for k in range(0,len(xnew)):  #mutate the parameters according to the bounds
                    if np.random.random() < Mut_rate:
                       # xnew[k] = (Bounds[1,k] - Bounds[0,k])*np.random.random()+Bounds[0,k]  
                        
                        
                        
                        xnew[k] = (Bounds[1,k]- Bounds[0,k]) + (10**(-4*np.random.random()))*np.random.random() + Bounds[0,k]
                '''
                changein = np.unique(np.random.randint(8, size=(Mut_rate*10)))
                
                
               
                for k in range(0,len(changein)):
                        tot_param[0][changein[k]] +=1
                        ind = changein[k]
                        xnew[ind] = x[ind] + np.random.normal(0,waver_param[0][ind])
                        if xnew[ind] > Bounds[1,ind]:
                            xnew[ind] = Bounds[1,ind] - abs(np.random.normal(0,.1))
                        if xnew[ind] < Bounds[0,ind]:
                            xnew[ind] = Bounds[0,ind] + abs(np.random.normal(0,.1))
                            
                fnew = f(xnew)  #evaluate function
                if fnew < fx or np.random.random(1) <= np.exp(fx - fnew): #test new function to keep
                    x = np.copy(xnew)
                    fx = np.copy(fnew)             #record this parameter combination
                    if fnew < fbest:        # update fbest and xbest
                        fbest = np.copy(fnew)
                        xbest = np.copy(xnew)
                        
                    for k in range(0,len(changein)): 
                        accept_param[0][changein[k]] += 1
                        
                    
                    
                if self.methastecounter >= 0:             
                    if j == N_Thin-1:
                        funchain = np.append(funchain,fx)
                        parchain = np.vstack([parchain, x])     
                        
            self.methasteprog.step(amount=1)  #update prog bar


            for n in range(len(accept_param[0])):
                
                acceptrate = accept_param[0][n]/tot_param[0][n]
                if tot_param[0][n] !=0:
                    
                    if acceptrate < .1:
                        waver_param[0][n] = waver_param[0][n]*1.1
                    if acceptrate > .4 :
                        waver_param[0][n] = waver_param[0][n]*.7
                        
                    if waver_param[0][n] > waver_max[0][n]:
                        waver_param[0][n] = waver_max[0][n]
                    if waver_param[0][n] < waver_min[0][n]:
                        waver_param[0][n] = waver_min[0][n]
            
         
            
            nowold = now
            now = int(np.floor(time.time()-starttime))
                
            if self.methastecounter < 1000:

                
                if nowold == now:
                    evalpersec= evalpersec + 1     
                    
                else:
                   # if len(predicttime) < 40:

                        
                    predicttime.append(int((N_Chain+N_Burn)/evalpersec))
                    predvar= (sum(predicttime[1:]))/(len(predicttime))                    
                    evalpersec = 1
                    
                
            m, s = divmod(now, 60)
            h, m = divmod(m, 60)
            curr = "%d:%02d:%02d" % (h, m, s)
            self.prog_clock.configure(text=curr)
            
            m, s = divmod(predvar, 60)
            h, m = divmod(m, 60)
            pred = "%d:%02d:%02d" % (h, m, s)            
            self.pred_clock.configure(text=pred)
            self.methasteprog_window.update()      

        self.methasteprog_window.destroy()        #destroy the GUI
        self.focus_force() #reurn focus to GUI
        return funchain, parchain, fbest, xbest     

    def pause_MH(self):
        if self.methaste_go == True:
            self.methaste_go == False
            self.pause_mh.config(text="Resume")
        if self.methaste_go == False:
            self.methaste_go == True
            self.pause_mh.config(text="Pause")
            
    def MetroplisHastings_FSP(self,Function,Init_Parameters,Bounds,N_Burn,N_Chain,N_Thin,Mut_rate,N):
        self.methasteprog_window = Toplevel()
        self.methasteprog_window.iconbitmap('qbio_icon_9nD_icon.ico')
        self.methasteprog_window.title('Met-Haste Progress')
        #self.methasteprog_window.focus_force()
        
        print Bounds
        methasteprog_frame = Frame(self.methasteprog_window)
        self.methasteprog = Progressbar(methasteprog_frame,orient=HORIZONTAL, length=320, mode='determinate')
        methasteprog_frame.pack()
        self.methasteprog_window.resizable(width=False,height=False)
        starttime = time.time()
        now = starttime
        evalpersec = 1
        self.methaste_go = True
        predvar = 0
        predicttime = []
        self.methasteprog.grid(row=0,column=0,columnspan=5)
        self.prog_clock = Label(methasteprog_frame,text='')
        self.prog_clock.grid(row=1,column=1,pady=2,padx=2,sticky=W)
        prog_label = Label(methasteprog_frame,text='Run Time:')
        pred_label = Label(methasteprog_frame,text='Prediction:')
        prog_label.grid(row=1,column=0,pady=2,padx=2,sticky=E)
        pred_label.grid(row=1,column=2,pady=2,padx=2,sticky=E)
        
        self.pause_mh = Button(methasteprog_frame,text='Pause',command = self.pause_MH)
        self.pause_mh.grid(row=2,column=1,columnspan=3,sticky=W)
        
        self.pred_clock = Label(methasteprog_frame,text='')
        self.pred_clock.grid(row=1,column=3,pady=2,padx=2,sticky=W)
        self.eval_amount = Label(methasteprog_frame,text='0/'+str(N_Burn+N_Chain))
        self.eval_amount.grid(row=1,column=4,pady=2,padx=2,sticky=W)       
        self.methasteprog["maximum"]=(N_Chain+N_Burn)*N_Thin
        self.methasteprog.grab_set()
        self.methastecounter = 0
        x0 = np.copy(Init_Parameters)
        xi = np.copy(x0)#self.make_binaryvec(x0,Bounds[0,:],Bounds[1,:],N,M) #convert the intial parameter guess to binary

        xbest = np.copy(xi) #intialize the xbest storage variable
       
        f = Function #define the function
        #xfloat = self.make_floatvec(x,Bounds[0,:],Bounds[1,:],N,M) #intialize xfloat
        fx,N = f(xi,N) #intialize fx
        fbest = np.copy(fx) #initalize fbest
        funchain = np.array([]) #intialize output variables
        parchain = Init_Parameters
        Ncounter = 0
        Nlist = np.array([]) 
        xi = Init_Parameters
        waver_param = np.ones((1,8))
        accept_param = np.zeros((1,8))
        tot_param = np.zeros((1,8))
        waver_max = np.ones((1,8))
        for i in range(0,len(waver_max[0])):
            waver_max[0][i] = np.sqrt(abs(Bounds[1,i]-Bounds[0,i]))
        
        waver_min = np.ones((1,8))
        for i in range(0,len(waver_max[0])):
            waver_min[0][i] = abs(Bounds[1,i]-Bounds[0,i])*.0005
        
        while self.methaste_go == True:

            for i in range(-N_Burn,N_Chain): #for the length of Nchain, record
                
                
                pre = time.time()
                for j in range(0,N_Thin):
                    
                    changein = np.unique(np.random.randint(8, size=(Mut_rate*10)))
                    
                    
                    xnew = np.copy(xi) #copy the binary vector
                    for k in range(0,len(changein)):
                            tot_param[0][changein[k]] +=1
                            ind = changein[k]
                            xnew[ind] = xi[ind] + np.random.normal(0,waver_param[0][ind])
                            if xnew[ind] > Bounds[1,ind]:
                                xnew[ind] = Bounds[1,ind] - abs(np.random.normal(0,.1))
                            if xnew[ind] < Bounds[0,ind]:
                                xnew[ind] = Bounds[0,ind] + abs(np.random.normal(0,.1))
                                    
    
                    #gaussian distrubted around each found parameter
                    #each have their own weight, each time you select index
                    #tuning parameter for each indices. 
                    #probablity (accept | xold) probability (xnew|xold)*P(accept|xold,xnew) dxold
                           
                    Nold = N
                    
                    fnew,N = f(xnew,N)
                    if N == Nold:
                        Ncounter += 1
                        if Ncounter > 4:
                            N = int(np.ceil(.7*N))
                
                    if fnew < fx or np.random.random(1) <= np.exp(fx - fnew):
                        
                        xi = np.copy(xnew)
                        fx = np.copy(fnew)
                        
                        for k in range(0,len(changein)): 
                            accept_param[0][changein[k]] += 1
                            
                        print accept_param, "update"
                            
                        
                        if fnew < fbest:
                            fbest = np.copy(fnew)
                            xbest = np.copy(xnew)
                          
                    if i >= 0:
                        
                        if j == N_Thin-1:
                        
                            funchain = np.append(funchain,fx)
                            parchain = np.vstack([parchain, xi])
                            Nlist = np.append(Nlist,N)
                            
                    self.methasteprog.step(amount=1)  #update prog bar
                    self.methasteprog_window.update()  
                if i >=0:
                    for n in range(len(accept_param[0])):
                        
                        acceptrate = accept_param[0][n]/tot_param[0][n]
                        if tot_param[0][n] !=0:
                           
                            if acceptrate < .1:
                                waver_param[0][n] = waver_param[0][n]*1.1
                            if acceptrate > .4 :
                                waver_param[0][n] = waver_param[0][n]*.7
                                
                            if waver_param[0][n] > waver_max[0][n]:
                                waver_param[0][n] = waver_max[0][n]
                            if waver_param[0][n] < waver_min[0][n]:
                                waver_param[0][n] = waver_min[0][n]
                            
        
               # tot_param = np.zeros((1,8))
                #accept_param = np.zeros((1,8))
                self.methastecounter += 1
                
                nowold = now
                
                now = int(np.floor(time.time()-starttime))
                evaltime = time.time() - pre
               
                predtime = (N_Chain+N_Burn)*evaltime
                
                '''
                if self.methastecounter < 1000:
    
    
                    if nowold == now:
                        evalpersec= evalpersec + 1     
                        
                    else:
                       # if len(predicttime) < 40:
    
                            
                        predicttime.append(int((N_Chain+N_Burn)/evalpersec))
                        predvar= (sum(predicttime[1:]))/(len(predicttime))                    
                        evalpersec = 1
                '''
                m, s = divmod(now, 60)
                h, m = divmod(m, 60)
                curr = "%d:%02d:%02d" % (h, m, s)
                self.prog_clock.configure(text=curr)
                
                m, s = divmod(predtime, 60)
                h, m = divmod(m, 60)
                pred = "%d:%02d:%02d" % (h, m, s)            
                self.pred_clock.configure(text=pred)
               
                self.eval_amount.configure(text=str(self.methastecounter)+'/'+str(N_Burn+N_Chain))
                
                if self.methastecounter == N_Burn+N_Chain:
                    self.methaste_go = False

        self.methasteprog_window.destroy()        #destroy the GUI
        self.focus_force() #reurn focus to GUI

                        
        return funchain, parchain, fbest, xbest, Nlist          
    ########################################################################
    # Simulated Annealing
    # weighted genetic algorithm search
    ########################################################################

    def SimulatedAnnealing(self,Function,Init_Parameters,Bounds,mut_rate,itermax):
        #create a focus force progress bar
        self.methasteprog_window = Toplevel()
        self.methasteprog_window.iconbitmap('qbio_icon_9nD_icon.ico')
        self.methasteprog_window.title('SA Progress')
        self.methasteprog_window.focus_force()
        methasteprog_frame = Frame(self.methasteprog_window)
        self.methasteprog = Progressbar(methasteprog_frame,orient=HORIZONTAL, length=241, mode='determinate')
        methasteprog_frame.pack()
        self.methasteprog_window.resizable(width=False,height=False)
        starttime = time.time()
        now = starttime
        evalpersec = 1
        evalmax = 1
        predvar = 0
        predicttime = []
        self.methasteprog.grid(row=0,column=0,columnspan=4)
        self.prog_clock = Label(methasteprog_frame,text='')
        self.prog_clock.grid(row=1,column=1,pady=2,padx=2,sticky=W)
        prog_label = Label(methasteprog_frame,text='Run Time:')
        pred_label = Label(methasteprog_frame,text='Prediction:')
        prog_label.grid(row=1,column=0,pady=2,padx=2,sticky=E)
        pred_label.grid(row=1,column=2,pady=2,padx=2,sticky=E)
        
        self.pred_clock = Label(methasteprog_frame,text='')
        self.pred_clock.grid(row=1,column=3,pady=2,padx=2,sticky=W)
        
        self.methasteprog.grid(row=0,column=0)
        self.methasteprog["maximum"]=itermax
        self.methasteprog.grab_set()
        
        #initialize parameters
        x0 = np.copy(Init_Parameters) #intial state
        funchain = np.array([])
        parchain = Init_Parameters    #intial param
        x = np.copy(x0)             #preallocate other vars
        xfloat = np.copy(x)
        xbest = np.copy(x)
        f = Function
        
        it = 0 #set up iteration
        fx = f(xfloat)
        fbest = np.copy(fx)   
        T=100  #temperature inital
        Tmin = 1e-6 #min Temp
        dT = .01 #Temp Change per step
        
        while T>Tmin and it < itermax:
            it +=1
            xnew = np.copy(x)
            for k in range(0,len(xnew)): #while T> Tmin mutate bases
                if np.random.random() < mut_rate:
                    xnew[k] = (Bounds[1,k] - Bounds[0,k])*np.random.random()+Bounds[0,k]              
                    
            fnew = f(xnew)
            
            if fnew <=fx or np.log(np.random.random(1)) <= ((fx-fnew)/T): #choose to keep the guesses based on temperature
                
                x = np.copy(xnew)
                fx = np.copy(fnew)
                if fnew < fbest:   #record and update fbest xbest
                    fbest = np.copy(fnew)
                    xbest = np.copy(xnew)
                           
                funchain = np.append(funchain,fx)
                parchain = np.vstack([parchain, xfloat])        
            T = T*(1-dT) #update temperature
            self.methasteprog.step(amount=1) #update progress bar
            nowold = now
            now = int(np.floor(time.time()-starttime))
            
            if self.methasteprog.getint() < 1000:
                if nowold == now:
                    evalpersec= evalpersec + 1     
    
                else:
                   # if len(predicttime) < 40:
    
                    predicttime.append(int((itermax)/evalpersec))
                    predvar= (sum(predicttime[1:]))/(len(predicttime))
                    evalpersec = 1
                                    
            m, s = divmod(now, 60)
            h, m = divmod(m, 60)
            curr = "%d:%02d:%02d" % (h, m, s)
            self.prog_clock.configure(text=curr)
            
            m, s = divmod(predvar, 60)
            h, m = divmod(m, 60)
            pred = "%d:%02d:%02d" % (h, m, s)            
            self.pred_clock.configure(text=pred)
            self.methasteprog_window.update()      

            self.methasteprog_window.update()      
        self.methasteprog_window.destroy()    #destroy progbar and return focus
        self.focus_force()
        floatbest = xbest        
        return floatbest,fbest,funchain,parchain
        
    
    def Second_Moment_Call(self):
        
        k12 = self.k12_entry.get()
        k23 = self.k23_entry.get()
        k21 = self.k21_entry.get()
        k32 = self.k32_entry.get()
        kr2 = self.kr2_entry.get()
        kr3 = self.kr3_entry.get()
        gamma =  self.gamma_entry.get()
        beta = self.beta_entry.get()
        
        parameters = np.array([float(k12),   float(k23),   float(k21),   float(k32), float(kr2), float(kr3),  float(gamma), float(beta)])
        
        #build time
        t = np.linspace(float(self.x0_t0.get()),float(self.x0_tf.get()),100*(float(self.x0_tf.get())-float(self.x0_t0.get())))
        case = self.case.get()
        x0 = [float(self.x0_S1.get()), float(self.x0_S2.get()), float(self.x0_S3.get()),float(self.x0_R.get()),0.,0.,0.,0.,0.,0.,0.,0.,0.,0.] #build x0
        
        if case == 0:        
                       #solve the solution for case 0
            checkfun=0                                                      #no need to check the function we dont use it
        else:
            #If the case is time variant check the user inputted function 
            checkfun = self.test_input_function(self.Input_entry.get())
            TimeVar = self.GeneModel.TimeVar_Create(self.Input_entry.get())
            #if the function given is ok, calculate the solution
        if checkfun == 0:
            if case == 0:
                solution = self.GeneModel.Solve_Moments(x0,t,parameters)  
            else:
                solution = self.GeneModel.Solve_Moments_TimeVar(x0,t,parameters,case,TimeVar)
            labels = ['<1>','<2>','<3>','<R>',	'sig 1','SIG_12','SIG_13','SIG_14','sig 2','SIG_23','SIG_24','sig 3','SIG_34',	'sig 4']
            self.fig.clf()
            self.subplot = self.fig.add_subplot(211)
            self.subplot.set_prop_cycle('color',self.colormap)
            
            for i in range(0,14):
                if i in [0,1,2,3,4,8,11,13]:
                    self.subplot.plot(t,solution[:,i],label=labels[i])

            self.subplot.legend(loc='best',fontsize=7)
            
            self.subplot = self.fig.add_subplot(212)
            self.subplot.set_prop_cycle('color',self.colormap)
            minvec = np.zeros((1,len(solution[:,3])))
            minvec = minvec.flatten()
            
            stdvec = np.sqrt(solution[:,13])
            for i in range(0,len(solution[:,3])):
                minvec[i] = max(solution[i,3]-solution[i,13],0)
            self.subplot.plot(t,solution[:,3],color=self.colormap[3])
            self.subplot.plot(t,solution[:,3]+solution[:,13],color=self.colormap[3],ls='dotted')
            self.subplot.plot(t,solution[:,3]+stdvec,color=self.colormap[3],ls='dashed')
            self.subplot.plot(t,minvec,color=self.colormap[3],ls='dashed') 
            
            ax2 = self.fig.gca()

            ax2.fill_between(t,minvec,solution[:,3]+solution[:,13],color='#e0e0f0')
            
            
            self.subplot.set_xlabel('$Time$')
            self.subplot.set_ylabel('$RNA$')
            self.figurecanvas.show()    
            self.fig.canvas.draw()


    ########################################################################
    # Prediction button callback
    # pulls entries and will run the required prediction
    # after it will update the plot
    ########################################################################
       
    def run_prediction(self):
        
        prediction_val = self.predict.get() #get the prediction type
        #get all parameters from entry boxes
        try:
            k12 = self.k12_entry.get()
            k23 = self.k23_entry.get()
            k21 = self.k21_entry.get()
            k32 = self.k32_entry.get()
            kr2 = self.kr2_entry.get()
            kr3 = self.kr3_entry.get()
            gamma =  self.gamma_entry.get()
            beta = self.beta_entry.get()
        
        #build parameters
            parameters = np.array([float(k12),   float(k23),   float(k21),   float(k32), float(kr2), float(kr3),  float(gamma), float(beta)])
        
        #build time
            t = np.linspace(float(self.x0_t0.get()),float(self.x0_tf.get()),10*(float(self.x0_tf.get())-float(self.x0_t0.get())))
        except:
            tkMessageBox.showerror(title = "Input Error", message = "Some of the parameter entries are not numbers or are empty, double check the parameters please.")
            return
        case = self.case.get()
        #Prediction val = 0 for ODE
        if prediction_val == 0:
            x0 = [float(self.x0_S1.get()), float(self.x0_S2.get()), float(self.x0_S3.get()),float(self.x0_R.get())] #build x0
            if case == 0:        
                solution = self.GeneModel.Gene_Dot(x0,t,parameters)             #solve the solution for case 0
                checkfun=0                                                      #no need to check the function we dont use it
            else:
                #If the case is time variant check the user inputted function 
                checkfun = self.test_input_function(self.Input_entry.get())
                TimeVar = self.GeneModel.TimeVar_Create(self.Input_entry.get())
                #if the function given is ok, calculate the solution
                if checkfun == 0:
                    solution = self.GeneModel.Gene_Dot_TimeVar(x0,t,parameters,case,TimeVar)
           
            if checkfun ==0: #plot the solution if it was calculated
                if self.holdon.get() == 0 or self.currplot.get() == 2:
                    self.fig.clf()
                    self.subplot = self.fig.add_subplot(111)
                    self.subplot.set_prop_cycle('color',self.colormap)
                
                
                self.subplot.plot(t,solution[:,0],label='S1')
                self.subplot.plot(t,solution[:,1],label='S2')
                self.subplot.plot(t,solution[:,2],label='S3')
                self.subplot.plot(t,solution[:,3],label='R')
                
                
                if self.holdon.get() == 0 or self.currplot.get() == 2:                #if hold on is off plot legend
                    self.subplot.legend(loc='best')
                    
                self.subplot.set_xlabel('$Time$')
                self.figurecanvas.show()
                self.fig.canvas.draw()
                self.currplot.set(1)
            
        if prediction_val == 1: #Prediction Val = 1 for SSA
            x0 = np.array([[float(self.x0_S1.get())], [float(self.x0_S2.get())], [float(self.x0_S3.get())],[float(self.x0_R.get())]])
            if case != 0:            #Check the user given function if case != 0
                checkfun = self.test_input_function(self.Input_entry.get())
            else:
                checkfun = 0
                
            if checkfun ==0:
                #plot the solution
                
                TimeVar = self.GeneModel.TimeVar_Create(self.Input_entry.get()) #create timeVar
                solution = self.SSA_u(x0,t,parameters,int(self.SSA_entry.get()),case,TimeVar) #run the average SSA's
                
                if self.holdon.get() == 0 or self.currplot.get() == 2: #if hold on is off delete the current plot and update
                    self.fig.clf()
                    self.subplot = self.fig.add_subplot(111)
                    self.subplot.set_prop_cycle('color',self.colormap)
                    
                self.subplot.plot(t,solution[0,:],label='S1')
                self.subplot.plot(t,solution[1,:],label='S2')
                self.subplot.plot(t,solution[2,:],label='S3')
                self.subplot.plot(t,solution[3,:],label='R')
                
                if self.holdon.get() == 0 or self.currplot.get() == 2:
                    self.subplot.legend(loc='best')
                    
                self.subplot.set_xlabel('$Time$')
                self.figurecanvas.show()
                self.fig.canvas.draw()
                self.currplot.set(1)
       
        if prediction_val == 2: #Prediction Val = 1 for FSP
           # x0 = [float(self.x0_S1.get()), float(self.x0_S2.get()), float(self.x0_S3.get()),float(self.x0_R.get())]
        
            #initualize current states
            N = 20 + int(float(self.x0_R.get())) #default 10 states at 0 RNA
            x0 = np.zeros((1,3*N)).T
            x0[int(float(self.x0_R.get()))] = int(float(self.x0_S1.get()))
            x0[int(float(self.x0_R.get()))+1] = int(float(self.x0_S2.get()))
            x0[int(float(self.x0_R.get()))+2] = int(float(self.x0_S3.get()))
            if case != 0:            #Check the user given function if case != 0
                checkfun = self.test_input_function(self.Input_entry.get())
            else:
                checkfun = 0
                
            if checkfun ==0:       #if the user input function is fine run the FSP      
                TimeVar = self.GeneModel.TimeVar_Create(self.Input_entry.get())
                t = np.linspace(float(self.x0_t0.get()),float(self.x0_tf.get()),(float(self.x0_tf.get())-float(self.x0_t0.get())))

                
                try:
                    error = float(self.FSP_entry.get())
                except:
                    tkMessageBox.showerror(title='FSP error entry',message ='The FSP error tolerance is not a number, please input a tolerance floating point number.')
                    return
              #solve at different Times
                    
                solution,N = self.GeneModel.FSPCalc_ODE_prealloc(parameters,TimeVar,case,N,x0,t,error)
        
                N = len(solution[0])/3 
                x1,solution1 = self.GeneModel.FSP_Return_RNA_Prob(solution[len(t)/4],N)
                x2,solution2 = self.GeneModel.FSP_Return_RNA_Prob(solution[len(t)/2],N)
                x3,solution3 = self.GeneModel.FSP_Return_RNA_Prob(solution[3*len(t)/4],N)
                x4,solution4 = self.GeneModel.FSP_Return_RNA_Prob(solution[-1],N)
                

                

                self.fig.clf()
                
                
                #define the subplots and update for the 4 time points
                self.subplot = self.fig.add_subplot(321) 
                self.subplot.set_prop_cycle('color',self.colormap)
                self.subplot.bar(x1,solution1,align='center',color = self.colormap[0],lw=0.0,width=1)
                self.subplot.set_xticks(x1)
                self.subplot.set_xlim(left=0)                
                self.subplot.set_ylabel('$Probability$')
                
                self.subplot.set_title('$T$ ' + '$_{' +str(int((t[-1:] - t[0])/4))+'}$')
                
                self.subplot = self.fig.add_subplot(322)
                self.subplot.set_prop_cycle('color',self.colormap)
                self.subplot.bar(x2,solution2,align='center',color = self.colormap[0],lw=0.0,width=1)
                self.subplot.set_xticks(x2)
                self.subplot.set_xlim(left=0)    
                self.subplot.set_title('$T$ ' +'$_{'+ str(int(2*(t[-1:] - t[0])/4))+'}$')
                
                self.subplot = self.fig.add_subplot(323)
                self.subplot.set_prop_cycle('color',self.colormap)
                self.subplot.bar(x3,solution3,align='center',color = self.colormap[0],lw=0.0,width=1)
                self.subplot.set_xticks(x3)
                self.subplot.set_xlim(left=0)    
                self.subplot.set_ylabel('$Probability$')
                
                self.subplot.set_title('$T$ ' + '$_{'+str(int(3*(t[-1:] - t[0])/4))+'}$')
                
                self.subplot = self.fig.add_subplot(324)
                self.subplot.set_prop_cycle('color',self.colormap)
                self.subplot.bar(x4,solution4,align='center',color = self.colormap[0],lw=0.0,width=1)
                self.subplot.set_xticks(x4)
                self.subplot.set_xlim(left=0)                 
                
                self.subplot.set_title('$T$ ' + '$_{'+str(int(t.item(-1)))+'}$')
                
                self.subplot = self.fig.add_subplot(325)
                           
                for i in range(0,len(t)):
                    x,RNAprob = self.GeneModel.FSP_Return_RNA_Prob(solution[i],N)
                    self.subplot.bar(x,RNAprob,align='center',color=cm.jet(1.*i/len(t)),lw=0.0,width=1)
                    
                self.subplot.set_xlabel('# $RNA$')
                self.subplot.set_title('$All$' + ' ' + '$Time$')
                self.subplot.set_ylabel('$Probability$')
                self.subplot.set_xticks(x)
                self.subplot.set_xlim(left=0)
                xbar=x
                
                

                self.subplot = self.fig.add_subplot(326)
                for i in range(1,len(t)):
                    x,y = self.GeneModel.FSP_smooth_graph(solution[i],N)
                    x = abs(x)
                    y = abs(y)
                    plt.plot(x,y,color=cm.jet(1.*i/len(t)))         

                self.subplot.set_xlabel('# $RNA$')        
                self.subplot.set_title('$All$' + ' ' + '$Time$')
                self.subplot.set_xticks(xbar)
                self.subplot.set_xlim(left=0)   
                
                
                self.fig.tight_layout()
                self.fig.canvas.draw()
                
                self.currplot.set(2)
                self.figurecanvas.show()
                
                
        
    def FSP_fit_fun(self,Parameters,case,t,x_in,mRNA_Hist,TimeVar,tol,N):

        x0 = np.zeros((1,3*N)).T
        x0[0] = x_in[0]
        x0[1] = x_in[1]
        x0[2] = x_in[2]
        x0[3] = x_in[3]
        
        
        soln,N = self.GeneModel.FSPCalc_ODE_prealloc(Parameters,TimeVar,case,N,x0,t,tol)
        logsum = 0
    
        for i in range(0,len(mRNA_Hist[0])):
            nhist = np.histogram(mRNA_Hist.T[i],bins = np.linspace(min(mRNA_Hist.T[i]),max(mRNA_Hist.T[i]),(max(mRNA_Hist.T[i])-min(mRNA_Hist.T[i])+1)))
          
            b,RNAprob = self.GeneModel.FSP_Return_RNA_Prob(soln[i],N)
            for j in range(0,len(nhist[0])):
                if j > len(RNAprob)-1:
                    RNAprob = np.append(RNAprob,1e-10)
                if RNAprob[j] <=0:
                    RNAprob[j] = 1e-10       
                logsum = logsum + nhist[0][j]*np.log(RNAprob[j])
         
        return -logsum,N
              

        
    def FSP_fit_fun_scipymin(self,Parameters,case,t,x_in,mRNA_Hist,TimeVar,tol):
        N = self.N
        x0 = np.zeros((1,3*N)).T
        x0[0] = x_in[0]
        x0[1] = x_in[1]
        x0[2] = x_in[2]
        x0[3] = x_in[3]
        
        
        soln,self.N = self.GeneModel.FSPCalc_ODE_prealloc(Parameters,TimeVar,case,N,x0,t,tol)
        logsum = 0
    
        for i in range(0,len(mRNA_Hist[0])):
            nhist = np.histogram(mRNA_Hist.T[i],bins = np.linspace(min(mRNA_Hist.T[i]),max(mRNA_Hist.T[i]),(max(mRNA_Hist.T[i])-min(mRNA_Hist.T[i])+1)))
          
            b,RNAprob = self.GeneModel.FSP_Return_RNA_Prob(soln[i],N)
            for j in range(0,len(nhist[0])):
                if j > len(RNAprob)-1:
                    RNAprob = np.append(RNAprob,1e-10)
                if RNAprob[j] <=0:
                    RNAprob[j] = 1e-10       
                logsum = logsum + nhist[0][j]*np.log(RNAprob[j])
         
        return -logsum
              

                


                
    ########################################################################
    # Fitting button callback
    # pulls entries and will run the required fit
    # after it will update the plot
    ########################################################################
           
        
    def Run_Fit(self):
        fitnum = self.fittype.get()
        odeorfsp = self.odeorfspvar.get()
        
        #check to see if a file or CSV is loaded
        try:
            len(self.mRNA_Hist[0,:])
            checkfun = 0
        except:
            tkMessageBox.showerror(title="CSV error",message="There is a problem with the loaded file. Double check that it is the correct csv.")
            checkfun = 1
            
        #if there is a csv loaded intialize parameters
        if checkfun ==0:
            k12 = self.k12_entry.get()
            k23 = self.k23_entry.get()
            k21 = self.k21_entry.get()
            k32 = self.k32_entry.get()
            kr2 = self.kr2_entry.get()
            kr3 = self.kr3_entry.get()
            gamma =  self.gamma_entry.get()
            beta = self.beta_entry.get()
            
            parameters = np.array([float(k12),   float(k23),   float(k21),   float(k32), float(kr2), float(kr3),  float(gamma),float(beta)])
            t = np.linspace(float(self.x0_t0.get()),float(self.x0_tf.get()),len(self.mRNA_Hist[0,:]))
            case = self.case.get()
            x0 = [float(self.x0_S1.get()), float(self.x0_S2.get()), float(self.x0_S3.get()),float(self.x0_R.get())]
            Data_Set=self.mRNA_Hist
            meanData_set = np.zeros((1,len(Data_Set[0,:])))
            standevData_set = np.zeros((1,len(Data_Set[0,:])))
            
            #calculate the mean and std of data_set before searches for faster computing
            for i in range(0,len(Data_Set[0,:])):
                meanData_set[0,i] = np.average(Data_Set[:,i])
                standevData_set[0,i] = np.std(Data_Set[:,i])
                
                
            #if ODE fit and default scipy min selected
            if odeorfsp==0:
                if fitnum ==0:
                    checkfun = 0
                    if case !=0:
                        checkfun = self.test_input_function(self.Input_entry.get())
                    if checkfun ==0:
                        TimeVar = self.GeneModel.TimeVar_Create(self.Input_entry.get())
                        Bounds = self.bounds
                        Bounds = ((self.bounds[0][0],self.bounds[1][0]),(self.bounds[0][1],self.bounds[1][1]),(self.bounds[0][2],self.bounds[1][2]),(self.bounds[0][3],self.bounds[1][3]),(self.bounds[0][4],self.bounds[1][4]),(self.bounds[0][5],self.bounds[1][5]),(self.bounds[0][6],self.bounds[1][6]),(self.bounds[0][7],self.bounds[1][7]))
                        fminmethodstring = self.fminmethod.get().replace(" ","")
                        self.Best_Par_ODE = self.GeneModel.MinDifferenceODE(case,t,x0,Data_Set,parameters,TimeVar,standevData_set,meanData_set,Bounds,fminmethodstring)
                        self.MinError.config(state=NORMAL) 
                        self.MinError.delete(0.0,'end')
                        self.MinError.insert(0.0,str(self.Best_Par_ODE.fun)[:-len(str(self.Best_Par_ODE.fun))/3])
                        self.MinError.config(state=DISABLED)                         
                        self.Best_Par_ODE=self.Best_Par_ODE.x 
                        
                        if case !=0:
                            
                            solution = self.GeneModel.Gene_Dot_TimeVar(x0,t,self.Best_Par_ODE,case,TimeVar)
                        else: 
                            
                            solution = self.GeneModel.Gene_Dot(x0,t,self.Best_Par_ODE)
                            
                        #ignore hold on for this function as it will show data + fit
                        self.fig.clf()
                        self.subplot = self.fig.add_subplot(111)
                        self.subplot.set_prop_cycle('color',self.colormap)
                        #self.subplot.plot(t,Data_Set[0,:],'k.',label='Data')
                        self.Plot_Raw_Data()
                        
                        self.subplot.plot(t,solution[:,0],label='S1')
                        self.subplot.plot(t,solution[:,1],label='S2')
                        self.subplot.plot(t,solution[:,2],label='S3')
                        self.subplot.plot(t,solution[:,3],label='R')
                        self.subplot.legend(loc='best')
                        self.subplot.set_xlabel('$Time$')
                        
                        self.figurecanvas.show()
                        self.fig.canvas.draw()
                        self.currplot.set(1)

                        
                        
                        
                if fitnum ==2: #if Met haste is selected
                    checkfun = 0
                    if case !=0:
                        checkfun = self.test_input_function(self.Input_entry.get())    
                    if checkfun ==0:
                        #intialize met hast params
                        Bounds=self.bounds
                        TimeVar = self.GeneModel.TimeVar_Create(self.Input_entry.get())                       
                        x_in = x0
                        
                        NBurn = int(self.N_Burn.get())
                        NChain = int(self.N_Chain.get())
                        
                        NThin = int(self.N_Thin.get())
                        MutRate=float(self.Mut_Rate.get())
                        self.methastecounter = IntVar()
                        
                        #def min function
                        Function = lambda x: self.DifferenceODE(x,case,t,x_in,Data_Set,TimeVar,standevData_set,meanData_set)
                        self.Best_Par_MetHaste = self.MetroplisHastings(Function,parameters,Bounds,NBurn,NChain,NThin,MutRate)
                        self.MinError.config(state=NORMAL) 
                        self.MinError.delete(0.0,'end')
                        self.MinError.insert(0.0,str(self.Best_Par_MetHaste[2])[:-len(str(self.Best_Par_MetHaste[2]))/3])
                        self.MinError.config(state=DISABLED)            
                        # once again clear plot we want to see function iterations
                        self.fig.clf()
                        self.subplot = self.fig.add_subplot(111)
                        self.subplot.set_prop_cycle('color',self.colormap)
                        self.subplot.plot(np.linspace(0,float(self.N_Chain.get()),float(self.N_Chain.get())),self.Best_Par_MetHaste[0],label='Feval')
                        self.subplot.legend(loc='best')
                        self.subplot.set_xlabel('Function Iterations')
                        self.figurecanvas.show()
                        self.fig.canvas.draw()
                        self.currplot.set(2)       

                        
                if fitnum ==1:
                    checkfun = 0
                    if case !=0:
                        checkfun = self.test_input_function(self.Input_entry.get())    
                    if checkfun ==0:
                        TimeVar = self.GeneModel.TimeVar_Create(self.Input_entry.get())                       
                        x_in = x0
                        Bounds = self.bounds
                        try:
                            #self.progress["maximum"]=float(self.SA_limit.get())
                            checkfun = 0
                        except:
                            tkMessageBox.showerror(title="Iteration Limit Error",message="Please double check the input for the Simulated Annealing iteration limit")
                            checkfun=1
                        if checkfun==0:
                            Function = lambda x: self.GeneModel.DifferenceODE(x,case,t,x_in,Data_Set,TimeVar,standevData_set,meanData_set)
                            self.Best_Par_SA = self.SimulatedAnnealing(Function,parameters,Bounds,.1,float(self.SA_limit.get()))
                            self.MinError.config(state=NORMAL) 
                            self.MinError.delete(0.0,'end')
                            self.MinError.insert(0.0,str(self.Best_Par_SA[1])[:-len(str(self.Best_Par_SA[1]))/3])
                            self.MinError.config(state=DISABLED)          
                            self.fig.clf()
                            self.subplot = self.fig.add_subplot(111)
                            self.subplot.set_prop_cycle('color',self.colormap)
                            self.subplot.plot(np.linspace(0,len(self.Best_Par_SA[2]),len(self.Best_Par_SA[2])),self.Best_Par_SA[2],label='Feval')
                            self.subplot.legend(loc='best')
                            self.subplot.set_xlabel('Function Iterations')
                            self.figurecanvas.show()  
                            self.fig.canvas.draw()
                            self.currplot.set(2)
                            
            if odeorfsp==1:
                
                if fitnum ==0:
                    checkfun = 0
                    if case !=0:
                        checkfun = self.test_input_function(self.Input_entry.get())
                    if checkfun ==0:
                        TimeVar = self.GeneModel.TimeVar_Create(self.Input_entry.get())
                        Bounds = self.bounds
                        Bounds = ((self.bounds[0][0],self.bounds[1][0]),(self.bounds[0][1],self.bounds[1][1]),(self.bounds[0][2],self.bounds[1][2]),(self.bounds[0][3],self.bounds[1][3]),(self.bounds[0][4],self.bounds[1][4]),(self.bounds[0][5],self.bounds[1][5]),(self.bounds[0][6],self.bounds[1][6]),(self.bounds[0][7],self.bounds[1][7]))
                        
                        
                        fminmethodstring = self.fminmethod.get().replace(" ","")
                        self.N = 20
                        x_in = x0
                        try:
                            tol = float(self.FittingFSP_entry.get())
                        except:
                            tkMessageBox.showerror(title='FSP error entry',message ='The FSP error tolerance is not a number, please input a tolerance floating point number.')
                        #def min function
                            return
                        Function = lambda x: self.FSP_fit_fun_scipymin(x,case,t,x_in,Data_Set,TimeVar,tol)

                        
                        self.Best_Par_ODE = self.GeneModel.MinDifferenceFSP(case,Bounds,fminmethodstring,Function,parameters)
 


                        self.MinError.config(state=NORMAL) 
                        self.MinError.delete(0.0,'end')
                        self.MinError.insert(0.0,str(self.Best_Par_ODE.fun)[:-len(str(self.Best_Par_ODE.fun))/3])
                        self.MinError.config(state=DISABLED)                         
                        self.Best_Par_ODE=self.Best_Par_ODE.x 
                        
                        if case !=0:
                            
                            solution = self.GeneModel.Gene_Dot_TimeVar(x0,t,self.Best_Par_ODE,case,TimeVar)
                        else: 
                            
                            solution = self.GeneModel.Gene_Dot(x0,t,self.Best_Par_ODE)
                            
                        #ignore hold on for this function as it will show data + fit
                        self.fig.clf()
                        self.subplot = self.fig.add_subplot(111)
                        self.subplot.set_prop_cycle('color',self.colormap)
                        #self.subplot.plot(t,Data_Set[0,:],'k.',label='Data')
                        self.Plot_Raw_Data()
                        
                        self.subplot.plot(t,solution[:,0],label='S1')
                        self.subplot.plot(t,solution[:,1],label='S2')
                        self.subplot.plot(t,solution[:,2],label='S3')
                        self.subplot.plot(t,solution[:,3],label='R')
                        self.subplot.legend(loc='best')
                        self.subplot.set_xlabel('$Time$')
                        
                        self.figurecanvas.show()
                        self.fig.canvas.draw()
                        self.currplot.set(1)
                
                
                
                
                
                if fitnum ==2:
                    checkfun = 0
                    if case !=0:
                        checkfun = self.test_input_function(self.Input_entry.get())    
                    if checkfun ==0:
                        #intialize met hast params
                        Bounds=self.bounds
                        TimeVar = self.GeneModel.TimeVar_Create(self.Input_entry.get())                       
                        x_in = x0
                        
                        
                        NBurn = int(self.N_Burn.get())
                        NChain = int(self.N_Chain.get())
                        
                        NThin = int(self.N_Thin.get())
                        MutRate=float(self.Mut_Rate.get())
                        self.methastecounter = IntVar()
  
                        try:
                            tol = float(self.FittingFSP_entry.get())
                        except:
                            tkMessageBox.showerror(title='FSP error entry',message ='The FSP error tolerance is not a number, please input a tolerance floating point number.')
                        #def min function
                            return
                       # pr.enable()
                        Function = lambda x,N: self.FSP_fit_fun(x,case,t,x_in,Data_Set,TimeVar,tol,N)
                        
                        self.Best_Par_MetHaste = self.MetroplisHastings_FSP(Function,parameters,Bounds,NBurn,NChain,NThin,MutRate,20)
                        #pr.disable()
                        #s = StringIO.StringIO()
                       # sortby = 'cumulative'
                       # ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
                       # ps.print_stats()
                       # print s.getvalue()
                        # once again clear plot we want to see function iterations
                        self.fig.clf()
                        self.subplot = self.fig.add_subplot(211)
                        self.subplot.set_prop_cycle('color',self.colormap)
                        self.subplot.plot(np.linspace(0,float(self.N_Chain.get()),float(self.N_Chain.get())),self.Best_Par_MetHaste[0],label='Feval')
                        self.subplot.legend(loc='best')
                        self.subplot.set_xlabel('Function Iterations')
  
                        
                        self.subplot = self.fig.add_subplot(212)
                        self.subplot.set_prop_cycle('color',self.colormap)
                        self.subplot.plot(np.linspace(0,float(self.N_Chain.get()),float(self.N_Chain.get())),self.Best_Par_MetHaste[4],label='FSP size')
                        self.subplot.legend(loc='best')
                        self.subplot.set_xlabel('Function Iterations')
                        self.subplot.set_ylabel('States')
                        self.figurecanvas.show()
                        self.fig.canvas.draw()
                        self.currplot.set(2)     
###############################################################################
# GUI running functions and backends
# all backend or GUI related fun defined after mathmatical ones
###############################################################################

    ########################################################################
    # Export the current plots to a text file
    # checks the plot type before exporting
    # Needs to be generalized
    ########################################################################       
    def Export_PlotVar(self):
        currentplottype = self.currplot.get() #get current plot type
        if currentplottype ==0:
            tkMessageBox.showerror(title="Plot Export Error",message="There is no current plot, please plot something before exporting")
            return
        if currentplottype ==1:
            ax = self.fig.gca()
            
            f = tkFileDialog.asksaveasfile(mode='w',defaultextension='.txt')
            for i in range(0,len(ax.lines)):
                line = ax.lines[i]
                if i ==0:
                    a = np.asarray(line.get_xdata())
                    b = np.asarray(line.get_ydata())      
                    np.savetxt(f,(a,b),delimiter=',',fmt='%1.4f',newline='\r\n')
                else:
                    b = np.asarray(line.get_ydata())
                    np.savetxt(f,(b),delimiter=',',fmt='%1.4f',newline='\r\n')
                    
        if currentplottype == 2:
            x=1



    def Change_cmap(self):
        
        currentplottype = self.currplot.get()
        if currentplottype ==3:
            tkMessageBox.showerror(title="Error: Cannot Change Colormap",message="There is no current plot, please plot something before exporting")
            return
        if currentplottype ==1 or currentplottype==0:
            # Create the popup window
            self.colormap_window = Toplevel(self)
            self.colormap_window.iconbitmap('qbio_icon_9nD_icon.ico')
            self.colormap_window.resizable(width=False, height=False)
            self.colormap_window.focus_force()
            self.colormap_frame = Frame(self.colormap_window,relief=RAISED)
            self.colormap_window.title("Change Plot Color cycle")
            self.colormap_window.grab_set()
            
            self.colormap_frame.pack()
           
            
            self.colormap1_select =  Radiobutton(self.colormap_frame,var = self.colormap_select_var, value = 0)           
            self.colormap1_select.grid(column=0,row=1,padx=2,pady=2)

            self.colormap2_select =  Radiobutton(self.colormap_frame,var = self.colormap_select_var, value = 1)           
            self.colormap2_select.grid(column=0,row=2,padx=2,pady=2)
            
            self.colormap3_select =  Radiobutton(self.colormap_frame,var = self.colormap_select_var, value = 2)           
            self.colormap3_select.grid(column=0,row=3,padx=2,pady=2)
            
            self.colormap4_select =  Radiobutton(self.colormap_frame,var = self.colormap_select_var, value = 3)           
            self.colormap4_select.grid(column=0,row=4,padx=2,pady=2)
            
            self.colormap5_select =  Radiobutton(self.colormap_frame,var = self.colormap_select_var, value = 4)           
            self.colormap5_select.grid(column=0,row=6,padx=2,pady=2)
            
            colormap1_1 = Canvas(self.colormap_frame,width=320,height=50)            
            colormap1_1.grid(column=1,row=1,padx=2,pady=2,columnspan =9)
            
            #create 4 canvases displaying the default colors
            self.colormap1_1_items = ['blue','red','green','magenta','cyan','yellow','black','#00FF00']
            for i in range(0,len( self.colormap1_1_items)):               
                colormap1_1.create_rectangle(i*40, 0, (i+1)*40, 40, fill=self.colormap1_1_items[i])
                
            colormap1_2 = Canvas(self.colormap_frame,width=320,height=50)            
            colormap1_2.grid(column=1,row=2,padx=2,pady=2,columnspan =9)
            self.colormap1_2_items = ['#000000','#888888','#444444','#CCCCCC','#222222','#AAAAAA','#666666','#EEEEEE']
            for i in range(0,len( self.colormap1_2_items)):               
                colormap1_2.create_rectangle(i*40, 0, (i+1)*40, 40, fill=self.colormap1_2_items[i])

            colormap1_3 = Canvas(self.colormap_frame,width=320,height=50)            
            colormap1_3.grid(column=1,row=3,padx=2,pady=2,columnspan =9)
            self.colormap1_3_items = ['#440154','#9FDA39','#365B8C','#FDE724','#4AC16D','#277E8E','#46327E','#9FD939']
            for i in range(0,len( self.colormap1_3_items)):               
                colormap1_3.create_rectangle(i*40, 0, (i+1)*40, 40, fill=self.colormap1_3_items[i])

            colormap1_4 = Canvas(self.colormap_frame,width=320,height=50)            
            colormap1_4.grid(column=1,row=4,padx=2,pady=2,columnspan =9)
            self.colormap1_4_items = ['#0c0787','#f48849','#8b0985','#f0f921','#5301a3','#dc5c64','#b93289','#febd2a']
            for i in range(0,len( self.colormap1_4_items)):               
                colormap1_4.create_rectangle(i*40, 0, (i+1)*40, 40, fill=self.colormap1_4_items[i])       
            
            #create the custom entry and have the current plot colors autofill in
            custom_label = Label(self.colormap_frame,text ="Custom, click to change colors")
            custom_label.grid(column=1,row=5,padx=2,pady=2,sticky=W,columnspan=5)
            self.colormap_custom_button_1 = Button(self.colormap_frame,width=4,height=2,bg =self.customcolormap[0],command=self.colorbutton_1 )
            self.colormap_custom_button_1.grid(column=1,row=6,padx=0,pady=2,sticky=W)

            self.colormap_custom_button_2 = Button(self.colormap_frame,width=4,height=2,bg = self.customcolormap[1],command=self.colorbutton_2 )
            self.colormap_custom_button_2.grid(column=2,row=6,padx=0,pady=2,sticky=W)
            self.colormap_custom_button_3 = Button(self.colormap_frame,width=4,height=2,bg = self.customcolormap[2],command=self.colorbutton_3 )
            self.colormap_custom_button_3.grid(column=3,row=6,padx=0,pady=2,sticky=W)
            self.colormap_custom_button_4 = Button(self.colormap_frame,width=4,height=2,bg = self.customcolormap[3],command=self.colorbutton_4 )
            self.colormap_custom_button_4.grid(column=4,row=6,padx=0,pady=2,sticky=W)
            self.colormap_custom_button_5 = Button(self.colormap_frame,width=4,height=2,bg = self.customcolormap[4],command=self.colorbutton_5 )
            self.colormap_custom_button_5.grid(column=5,row=6,padx=0,pady=2,sticky=W)
            
            self.colormap_custom_button_6 = Button(self.colormap_frame,width=4,height=2,bg = self.customcolormap[5],command=self.colorbutton_6 )
            self.colormap_custom_button_6.grid(column=6,row=6,padx=0,pady=2,sticky=W)
            self.colormap_custom_button_7 = Button(self.colormap_frame,width=4,height=2,bg =self.customcolormap[6],command=self.colorbutton_7 )
            self.colormap_custom_button_7.grid(column=7,row=6,padx=0,pady=2,sticky=W)
            self.colormap_custom_button_8 = Button(self.colormap_frame,width=4,height=2,bg = self.customcolormap[7],command=self.colorbutton_8 )
            self.colormap_custom_button_8.grid(column=8,row=6,padx=0,pady=2,sticky=W)
            self.set_cmap = Button(self.colormap_frame,width=15,text='Set Color Map', command = self.Update_Cmap)
            self.set_cmap.grid(column=5,row=7,padx=2,pady=4,columnspan=3)

    def colorbutton_1(self):
        c = askcolor()
        self.colormap_custom_button_1.config(bg=c[1])
        self.colormap_select_var.set(4)

    def colorbutton_2(self):
        c = askcolor()
        self.colormap_custom_button_2.config(bg=c[1])
        self.colormap_select_var.set(4)
    def colorbutton_3(self):
        c = askcolor()
        self.colormap_custom_button_3.config(bg=c[1])
        self.colormap_select_var.set(4)
    def colorbutton_4(self):
        c = askcolor()
        self.colormap_custom_button_4.config(bg=c[1])
        self.colormap_select_var.set(4)
    def colorbutton_5(self):
        c = askcolor()
        self.colormap_custom_button_5.config(bg=c[1])
        self.colormap_select_var.set(4)
    def colorbutton_6(self):
        c = askcolor()
        self.colormap_custom_button_6.config(bg=c[1])
        self.colormap_select_var.set(4)
    def colorbutton_7(self):
        c = askcolor()
        self.colormap_custom_button_7.config(bg=c[1])
        self.colormap_select_var.set(4)
    def colorbutton_8(self):
        c = askcolor()
        self.colormap_custom_button_8.config(bg=c[1])
        self.colormap_select_var.set(4)

    ########################################################################
    # Colormap Callback
    # checks and updates custom inputted colormap
    ########################################################################

    def Update_Cmap(self):      
        
        #if default no need to check the colormap update.
        if self.colormap_select_var.get() == 0:
            self.colormap = self.colormap1_1_items
            
            self.colormap_window.destroy()
        if self.colormap_select_var.get() == 1:
            self.colormap = self.colormap1_2_items
            self.colormap_window.destroy()
        if self.colormap_select_var.get() == 2:
            self.colormap = self.colormap1_3_items
            self.colormap_window.destroy()
        if self.colormap_select_var.get() == 3:
            self.colormap = self.colormap1_4_items
            self.colormap_window.destroy()
        if self.colormap_select_var.get() == 4:
            self.colormap = [self.colormap_custom_button_1.cget('bg'),self.colormap_custom_button_2.cget('bg'),self.colormap_custom_button_3.cget('bg'),self.colormap_custom_button_4.cget('bg'),self.colormap_custom_button_5.cget('bg'),self.colormap_custom_button_6.cget('bg'),self.colormap_custom_button_7.cget('bg'),self.colormap_custom_button_8.cget('bg')]
            self.customcolormap = [self.colormap_custom_button_1.cget('bg'),self.colormap_custom_button_2.cget('bg'),self.colormap_custom_button_3.cget('bg'),self.colormap_custom_button_4.cget('bg'),self.colormap_custom_button_5.cget('bg'),self.colormap_custom_button_6.cget('bg'),self.colormap_custom_button_7.cget('bg'),self.colormap_custom_button_8.cget('bg')]

        self.subplot.set_prop_cycle('color',self.colormap) #update cmap
        
        self.colormap_window.destroy()
        

    def openfile(self): #function to open the csv file with the mRNA_hist
           self.filename = askopenfilename(parent=self.OpenFile)
           if self.filename != "":
               try:
                   self.mRNA_Hist = self.GeneModel.Import_mRNA(self.filename)
               except:
                   tkMessageBox.showerror(title='Import Error',message='Something went wrong importing this file, please use a csv file or a data file with commas as the delimiter')
                   return
               self.FileEntry.config(state=NORMAL)
               self.FileEntry.delete(0.0,'end')
               self.FileEntry.insert(INSERT,self.filename+" loaded")
               self.FileEntry.config(state=DISABLED)

           
    ########################################################################
    # Define a function to create new figures
    ########################################################################
                             
    def Open_Fig(self):
        
#        self.plottingwindow = Toplevel(self)    #make a new window
#        self.plottingwindow.iconbitmap('qbio_icon_9nD_icon.ico')
#        self.plottingFrame=Frame(self.plottingwindow,relief=RAISED)
#        self.plottingFrame.pack(fill=BOTH)
#        self.plottingwindow.geometry("700x800+420+10")
#        self.plottingwindow.title('Figure')
        
      
        self.Fitting = IntVar()
        self.Fitting.set(0)
        #define fiture size and color
        
        self.fig = plt.figure()
        self.ax = self.fig.gca()
        self.ax.set_prop_cycle('color',self.colormap)
        #create call back to call figure plots
        self.subplot = self.fig.add_subplot(111)
        

        self.fig.set_facecolor('#DCDCDC') #imbed the canvase and define it as a tkfunction
        #self.figurecanvas = FigureCanvasTkAgg(self.fig,master=self.plottingFrame)
        self.fig.show()
        self.figurecanvas = self.fig
        #self.figurecanvas.get_tk_widget().pack(side=TOP,fill=BOTH,expand=1,padx=2,pady=2) #add default toolbar
        #toolbar = NavigationToolbar2TkAgg(self.figurecanvas, self.plottingFrame)
        #toolbar.update()
        #self.figurecanvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)        

    ########################################################################
    # Met Haste parameter space logic
    ########################################################################

    def par_space_setup(self,var1,used,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist):

        addon = 0 #logicl loop using plot flag var for each parameter check button
        while addon == 0:
            if k12plot ==1 and 0 not in used: 
                if addon == 0:
                    var1 = k12list
                    var1str = "k{12}"
                    used.append(0)
                    addon = 1
            if k23plot ==1 and 1 not in used:
                if addon == 0:
                    var1 = k23list
                    var1str = "k{23}"
                    used.append(1)
                    addon = 1
            if k21plot ==1 and 2 not in used:
                if addon == 0:
                    var1 = k21list
                    var1str = "k{21}"
                    used.append(2)
                    addon = 1
            if k32plot ==1 and 3 not in used: 
                if addon == 0:                
                    var1 = k32list
                    var1str = "k{32}"
                    used.append(3)
                    addon = 1
            if kr2plot ==1 and 4 not in used :
                if addon == 0:
                    var1 = kr2list
                    var1str = "k{r2}"
                    used.append(4)
                    addon = 1
            if kr3plot ==1 and 5 not in used:
                if addon == 0:
                    var1 = kr3list
                    var1str = "k{r3}"
                    used.append(5)
                    addon = 1
            if gmmplot ==1 and 6 not in used:
                if addon == 0:
                    var1 = gmmlist
                    var1str = "gamma"
                    used.append(6)
                    addon = 1
            if betplot ==1 and 7 not in used:
                if addon == 0:
                    var1 = betlist
                    var1str = "beta"
                    used.append(7)
                    addon = 1       
        #return the next variable and name not already used for plotting
        return var1,var1str,used


    def Plot_Raw_Data(self):
        
        try:
            len(self.mRNA_Hist[:])
            Data_Set = self.mRNA_Hist
            checkfun = 0
        except:
            tkMessageBox.showerror(title="CSV error",message="There is a problem with the loaded file. Double check that it is the correct csv.")
            checkfun = 1
            
        try:
            
            t = np.linspace(float(self.x0_t0.get()),float(self.x0_tf.get()),len(Data_Set[0,:]))
            
            checkfun=0
        except:
            tkMessageBox.showerror(title="Input error",message="There is a problem with the time inputs, double check the inputs for t0 and tf")
            checkfun = 1
            
        if checkfun ==0:
            
            
            meanData_set = np.zeros((1,len(Data_Set[0,:])))
            varData_set = np.zeros((1,len(Data_Set[0,:])))
            stdData_set = np.zeros((1,len(Data_Set[0,:])))
            for i in range(0,len(Data_Set[0,:])):
                meanData_set[0,i] = np.average(Data_Set[:,i])
                varData_set[0,i] = np.var(Data_Set[:,i])
                stdData_set[0,i] = np.std(Data_Set[:,i])
                
            self.fig.clf()
            self.subplot = self.fig.add_subplot(111)
            meanData_set = meanData_set.flatten()
            varData_set = varData_set.flatten()
            stdData_set = stdData_set.flatten()
            
            
            self.subplot.plot(t,meanData_set,'k.',label='$Data$')
            
            self.subplot.plot(t,meanData_set+stdData_set,'-',label='$- standev$',color='#c0c0c0')
            minvec = np.zeros(len(meanData_set))
            for i in range(0,len(meanData_set)):
                minvec[i] = max(meanData_set[i]-stdData_set[i],0)
            self.subplot.fill_between(t,minvec,meanData_set+stdData_set,color='#efefef')
            self.subplot.plot(t,minvec,'-',label='$+ standev$',color='#c0c0c0')
            self.subplot.legend(loc='best')


            self.figurecanvas.show()
            self.fig.canvas.draw()
            
    
    
    def Plot_Raw_Data_Hist(self):
        
        try:
            len(self.mRNA_Hist[:])
            
            checkfun = 0
        except:
            tkMessageBox.showerror(title="CSV error",message="There is a problem with the loaded file. Double check that it is the correct csv.")
            checkfun = 1
        if checkfun ==0:
            Data_Set = self.mRNA_Hist
            t = np.linspace(0,len(self.mRNA_Hist[1]),len(self.mRNA_Hist[1]))
            meanData_set = np.zeros((1,len(Data_Set[0,:])))
            varData_set = np.zeros((1,len(Data_Set[0,:])))
            stdData_set = np.zeros((1,len(Data_Set[0,:])))
            for i in range(0,len(Data_Set[0,:])):
                meanData_set[0,i] = np.average(Data_Set[:,i])
                varData_set[0,i] = np.var(Data_Set[:,i])
                stdData_set[0,i] = np.std(Data_Set[:,i])
                
            self.fig.clf()
            self.subplot = self.fig.add_subplot(211)
            meanData_set = meanData_set.flatten()
            varData_set = varData_set.flatten()
            stdData_set = stdData_set.flatten()
            
            
            self.subplot.plot(t,meanData_set,'k.',label='$Data$')
            self.subplot.plot(t,meanData_set-stdData_set,'-',label='$+ standev$',color='#c0c0c0')
            self.subplot.plot(t,meanData_set+stdData_set,'-',label='$- standev$',color='#c0c0c0')
            minvec = np.zeros(len(meanData_set))
            for i in range(0,len(meanData_set)):
                minvec[i] = max(meanData_set[i]-stdData_set[i],0)
            self.subplot.fill_between(t,minvec,meanData_set+stdData_set,color='#efefef')
                
            self.subplot.legend(loc='best')

                  
            self.subplot = self.fig.add_subplot(212)
            for i in range(0,len(t)):

                self.subplot.hist(self.mRNA_Hist.T[i],bins = np.linspace(min(self.mRNA_Hist.T[i]),max(self.mRNA_Hist.T[i]),(max(self.mRNA_Hist.T[i])-min(self.mRNA_Hist.T[i])+1)), color = cm.jet(1.*i/len(t)))
               # self.subplot.hist(nhist[1],nhist[0],align='center',color=cm.jet(1.*i/len(t)),lw=0.0,width=1)
                
            self.subplot.set_xlabel("# Data")
            self.fig.tight_layout()
            self.figurecanvas.show()
            self.fig.canvas.draw()
            
            
    def PlotParameterSpace(self):

        k12plot = self.k12checkvar.get()
        k21plot = self.k21checkvar.get()
        k23plot = self.k23checkvar.get()
        k32plot = self.k32checkvar.get()
        kr2plot = self.kr2checkvar.get()
        kr3plot = self.kr3checkvar.get()
        gmmplot = self.gmmcheckvar.get()
        betplot = self.betcheckvar.get()
        
        #check to see if a met haste search has been preformed
        try:
            len(self.Best_Par_MetHaste[1])
            checkfun = 0
        except:
            tkMessageBox.showerror(title="MetHaste Error",message="There is no Metropolis-Hasting fitting data currently, please Try Fit")
            checkfun =1
        #check to see if more than 1 checkbutton is selected
        if sum([k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot]) in [0,1] and checkfun==0:
            tkMessageBox.showerror(title="Parameter Error",message="There are not enough parameters checked to compare")
            checkfun = 1
            
        if checkfun==0:
            parchain = self.Best_Par_MetHaste[1]
            
            k12list = []
            k21list = []
            k23list = []
            k32list = []
            kr2list = []
            kr3list = []
            gmmlist = []
            betlist = []
            
            
            for i in range(0,len(parchain)):
                if k12plot ==1:
                    k12list.append(parchain[i][0])
                if k23plot ==1:
                    k23list.append(parchain[i][1])
                if k21plot ==1:
                    k21list.append(parchain[i][2])
                if k32plot ==1:
                    k32list.append(parchain[i][3])
                if kr2plot ==1:
                    kr2list.append(parchain[i][4])
                if kr3plot ==1:
                    kr3list.append(parchain[i][5])
                if gmmplot ==1:
                    gmmlist.append(parchain[i][6])
                if betplot ==1:
                    betlist.append(parchain[i][7])
  
            self.fig.clf() 
            self.currplot.set(2)
            self.ParSpace_Plot(self,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist)             
                                
                        
                                                
                        
                
        
    def activate_fittype_box(self): #function to gray out SSA and FSP boxes
        
        if self.fittype.get() == 0: #reconfigure each box 1=SSA selected
            self.MetHasteParSpace.configure(state=DISABLED)
            self.N_Burn.configure(state=DISABLED)
            self.N_Thin.configure(state=DISABLED)
            self.N_Chain.configure(state=DISABLED)
            self.Mut_Rate.configure(state=DISABLED)
            self.SA_limit.config(state=DISABLED)
            
            
            self.k12check.configure(state=DISABLED)
            self.k21check.configure(state=DISABLED)
            self.k23check.configure(state=DISABLED)
            self.k32check.configure(state=DISABLED)
            self.kr2check.configure(state=DISABLED)
            self.kr3check.configure(state=DISABLED)
            self.gmmcheck.configure(state=DISABLED)
            self.betcheck.configure(state=DISABLED)
            
        if self.fittype.get() == 1: #reconfigure each box 1=SSA selected
            self.MetHasteParSpace.configure(state=DISABLED)
            self.N_Burn.configure(state=DISABLED)
            self.N_Thin.configure(state=DISABLED)
            self.N_Chain.configure(state=DISABLED)
            self.Mut_Rate.configure(state=DISABLED)
            self.SA_limit.config(state=NORMAL)
            self.k12check.configure(state=DISABLED)
            self.k21check.configure(state=DISABLED)
            self.k23check.configure(state=DISABLED)
            self.k32check.configure(state=DISABLED)
            self.kr2check.configure(state=DISABLED)
            self.kr3check.configure(state=DISABLED)
            self.gmmcheck.configure(state=DISABLED)
            self.betcheck.configure(state=DISABLED)
            
        if self.fittype.get() == 2: # 2 FSP selected
            self.MetHasteParSpace.configure(state=NORMAL)
            self.N_Burn.configure(state=NORMAL)
            self.N_Thin.configure(state=NORMAL)
            self.N_Chain.configure(state=NORMAL)
            self.Mut_Rate.configure(state=NORMAL)
            self.SA_limit.config(state=DISABLED)
            self.k12check.configure(state=NORMAL)
            self.k21check.configure(state=NORMAL)
            self.k23check.configure(state=NORMAL)
            self.k32check.configure(state=NORMAL)
            self.kr2check.configure(state=NORMAL)
            self.kr3check.configure(state=NORMAL)
            self.gmmcheck.configure(state=NORMAL)
            self.betcheck.configure(state=NORMAL)
                
            
    def test_input_function(self,string):
        #Function to check the input function and display an error if it can not evaluate
        checkvar = 0
        
        try:
            TimeVar = self.GeneModel.TimeVar_Create(self.Input_entry.get())
            TimeVar(1)
        except:
            tkMessageBox.showerror(title="Input Error",message="There is a problem with the Input Function. Double check the function syntax.")
            
            checkvar = 1
        #Return a variable to continue the calculations
        return checkvar
        
            
    ########################################################################
    # change bounds 
    # opens a new top level window to let the user change parameter bounds
    ########################################################################
    def Change_Bounds(self):
        
        self.boundswindow = Toplevel(self) #toplevel opens a new window
        self.boundswindow.iconbitmap('qbio_icon_9nD_icon.ico') #set icon to qbio
        
        #dont let the user use the main window or resize the popup
        self.boundswindow.resizable(width=False, height=False) 
        self.boundswindow.focus_force()
        
        self.BoundsFrame=Frame(self.boundswindow,relief=RAISED)
        self.BoundsFrame.pack(fill=BOTH) #create a frame for later GUI elements
        
        #define the popup title and size
        #self.boundswindow.geometry("247x246+300+10")
        self.boundswindow.title('Change Bounds')
        #focus attention on this window
        self.boundswindow.grab_set()
        self.boundswindow.minsize(width=247,height=246)
        
        # make 3 columns and titled param, min, max
        boundslabel1 = Label(self.BoundsFrame,text="Parameters")
        boundslabel2 = Label(self.BoundsFrame,text="Min")
        boundslabel3 = Label(self.BoundsFrame,text="Max")
        boundslabel1.grid(row=0,column=0,padx=2,pady=2)
        boundslabel2.grid(row=0,column=1,padx=2,pady=2)
        boundslabel3.grid(row=0,column=2,padx=2,pady=2)
        
        #create entries for all parameter min maxes        
        self.k12_min = Entry(self.BoundsFrame,width=7,justify=RIGHT) 
        self.k12_max = Entry(self.BoundsFrame,width=7,justify=RIGHT) 
        self.k21_min = Entry(self.BoundsFrame,width=7,justify=RIGHT) 
        self.k21_max = Entry(self.BoundsFrame,width=7,justify=RIGHT) 
        self.k23_min = Entry(self.BoundsFrame,width=7,justify=RIGHT) 
        self.k23_max = Entry(self.BoundsFrame,width=7,justify=RIGHT)
        self.k32_min = Entry(self.BoundsFrame,width=7,justify=RIGHT) 
        self.k32_max = Entry(self.BoundsFrame,width=7,justify=RIGHT) 
        self.kr2_min = Entry(self.BoundsFrame,width=7,justify=RIGHT) 
        self.kr2_max = Entry(self.BoundsFrame,width=7,justify=RIGHT) 
        self.kr3_min = Entry(self.BoundsFrame,width=7,justify=RIGHT) 
        self.kr3_max = Entry(self.BoundsFrame,width=7,justify=RIGHT) 
        self.gamma_min = Entry(self.BoundsFrame,width=7,justify=RIGHT) 
        self.gamma_max = Entry(self.BoundsFrame,width=7,justify=RIGHT) 
        self.beta_min = Entry(self.BoundsFrame,width=7,justify=RIGHT) 
        self.beta_max = Entry(self.BoundsFrame,width=7,justify=RIGHT) 
        
        #label all the min-maxes        
        k12_label = Label(self.BoundsFrame,width=3,text="k" +u"\u2081"+u"\u2082") 
        k21_label = Label(self.BoundsFrame,width=3,text="k" +u"\u2082"+u"\u2081") 
        k23_label = Label(self.BoundsFrame,width=3,text="k" +u"\u2082"+u"\u2083") 
        k32_label = Label(self.BoundsFrame,width=3,text="k" +u"\u2083"+u"\u2082") 
        kr2_label = Label(self.BoundsFrame,width=3,text="k" + u"\u1D63"+ u"\u2082") 
        kr3_label = Label(self.BoundsFrame,width=3,text="k"+ u"\u1D63"+ u"\u2083") 
        gamma_label = Label(self.BoundsFrame,width=3,text=" γ ") 
        beta_label = Label(self.BoundsFrame,width=3,text=" β ")
        
        #pack all these elements 
        k12_label.grid(row=1,column=0,padx=2,pady=2)
        self.k12_min.grid(row=1,column=1,padx=2,pady=2)
        self.k12_max.grid(row=1,column=2,padx=2,pady=2)
        
        k21_label.grid(row=2,column=0,padx=2,pady=2)
        self.k21_min.grid(row=2,column=1,padx=2,pady=2)
        self.k21_max.grid(row=2,column=2,padx=2,pady=2) 
        
        k23_label.grid(row=3,column=0,padx=2,pady=2)
        self.k23_min.grid(row=3,column=1,padx=2,pady=2)
        self.k23_max.grid(row=3,column=2,padx=2,pady=2)
        
        k32_label.grid(row=4,column=0,padx=2,pady=2)
        self.k32_min.grid(row=4,column=1,padx=2,pady=2)
        self.k32_max.grid(row=4,column=2,padx=2,pady=2)

        k23_label.grid(row=5,column=0,padx=2,pady=2)
        self.k23_min.grid(row=5,column=1,padx=2,pady=2)
        self.k23_max.grid(row=5,column=2,padx=2,pady=2)
        
        kr2_label.grid(row=6,column=0,padx=2,pady=2)
        self.kr2_min.grid(row=6,column=1,padx=2,pady=2)
        self.kr2_max.grid(row=6,column=2,padx=2,pady=2)

        kr3_label.grid(row=7,column=0,padx=2,pady=2)
        self.kr3_min.grid(row=7,column=1,padx=2,pady=2)
        self.kr3_max.grid(row=7,column=2,padx=2,pady=2)
        
        gamma_label.grid(row=8,column=0,padx=2,pady=2)
        self.gamma_min.grid(row=8,column=1,padx=2,pady=2)
        self.gamma_max.grid(row=8,column=2,padx=2,pady=2)
        
        beta_label.grid(row=9,column=0,padx=2,pady=2)
        self.beta_min.grid(row=9,column=1,padx=2,pady=2)
        self.beta_max.grid(row=9,column=2,padx=2,pady=2)
        
       # space_label = Label(self.BoundsFrame,width=1,text=" ")
        #space_label.grid(row=10,column=0,pady=2,padx=1)
        
        #define a button to update all the bounds
        self.update_bounds = Button(self.BoundsFrame,text="Update Bounds",width=12,height=3,command=self.Update_Bounds)
        self.update_bounds.grid(row=11,column=1,rowspan=1,columnspan=3,padx=2,pady=2)
        
        #insert the current bounds into entry boxes for the user
        self.k12_min.insert(0,self.bounds[0,0])
        self.k21_min.insert(0,self.bounds[0,2])
        self.k23_min.insert(0,self.bounds[0,1])
        self.k32_min.insert(0,self.bounds[0,3])
        self.kr2_min.insert(0,self.bounds[0,4])
        self.kr3_min.insert(0,self.bounds[0,5])
        self.gamma_min.insert(0,self.bounds[0,6])
        self.beta_min.insert(0,self.bounds[0,7])
        
        self.k12_max.insert(0,self.bounds[1,0])
        self.k21_max.insert(0,self.bounds[1,2])
        self.k23_max.insert(0,self.bounds[1,1])
        self.k32_max.insert(0,self.bounds[1,3])
        self.kr2_max.insert(0,self.bounds[1,4])
        self.kr3_max.insert(0,self.bounds[1,5])
        self.gamma_max.insert(0,self.bounds[1,6])  
        self.beta_max.insert(0,self.bounds[1,7])
    
    def Update_Bounds(self):
        try: 
            #try to update bounds as floats
            self.bounds[0,0] = float(self.k12_min.get())
            self.bounds[0,2] = float(self.k21_min.get())
            self.bounds[0,1] = float(self.k23_min.get())
            self.bounds[0,3] = float(self.k32_min.get())
            self.bounds[0,4] = float(self.kr2_min.get())
            self.bounds[0,5] = float(self.kr3_min.get())
            self.bounds[0,6] = float(self.gamma_min.get())
            self.bounds[0,7] = float(self.beta_min.get())
            self.bounds[1,0] = float(self.k12_max.get())
            self.bounds[1,2] = float(self.k21_max.get())
            self.bounds[1,1] = float(self.k23_max.get())
            self.bounds[1,3] = float(self.k32_max.get())
            self.bounds[1,4] = float(self.kr2_max.get())
            self.bounds[1,5] = float(self.kr3_max.get())
            self.bounds[1,6] = float(self.gamma_max.get())
            self.bounds[1,7] = float(self.beta_max.get())
            
            #update the user with a message
            tkMessageBox.showinfo(message="Parameter Bounds Updated")
            #destroy the window and reset focus to the main GUI
            self.boundswindow.destroy()
            self.focus_force()
            
        except:
            #if failed error message
            tkMessageBox.showerror(title="Input Error",message="Something that is not a number was given as a bound, double check all the bounds")
        
    def Sim_SSA(self):
        self.ask_step = Toplevel()
        
        self.ask_step.iconbitmap('qbio_icon_9nD_icon.ico')
        self.ask_step.title('Step Size')
        
        self.timestep_label = Label(self.ask_step,text='How many time points?')
        self.timestep_label.grid(column=0,row=0,padx=2,pady=2)
        
        self.timestep = Entry(self.ask_step,width=13)
        self.timestep.grid(column=0,row=1,padx=2,pady=2)
        
        self.timestep.bind("<Return>",self.continue_sim)      

                
            
            
    def continue_sim(self,event):
        
        if self.timestep.get() !="":
            try:
                
                int(self.timestep.get())*2
            except:
                return
            ind = self.timestep.get()
            self.timestep.unbind("<Return>")
            
            self.ask_step.destroy()
        
            k12 = self.k12_entry.get()
            k23 = self.k23_entry.get()
            k21 = self.k21_entry.get()
            k32 = self.k32_entry.get()
            kr2 = self.kr2_entry.get()
            kr3 = self.kr3_entry.get()
            gamma =  self.gamma_entry.get()
            beta = self.beta_entry.get()
            
            #build parameters
            parameters = np.array([float(k12),   float(k23),   float(k21),   float(k32), float(kr2), float(kr3),  float(gamma), float(beta)])
            
            #build time
            t = np.linspace(float(self.x0_t0.get()),float(self.x0_tf.get()),(float(self.x0_tf.get())-float(self.x0_t0.get())))
            case = self.case.get()
            Nruns = int(self.SSA_entry.get())
            x0 = np.array([[float(self.x0_S1.get())], [float(self.x0_S2.get())], [float(self.x0_S3.get())],[float(self.x0_R.get())]])
            if case != 0:            #Check the user given function if case != 0
                checkfun = self.test_input_function(self.Input_entry.get())
            else:
                checkfun = 0
                
            if checkfun ==0:
                #plot the solution
                
                TimeVar = self.GeneModel.TimeVar_Create(self.Input_entry.get()) #create timeVar
    
            
                nt = (len(t),Nruns) #Preallocate the matrices
                S1 = np.zeros(nt)
                S2 = np.zeros(nt)
                S3 = np.zeros(nt)
                R = np.zeros(nt)    
    
            self.methasteprog_window = Toplevel()
            self.methasteprog_window.iconbitmap('qbio_icon_9nD_icon.ico')
            self.methasteprog_window.title('SSA Progress')
        #self.methasteprog_window.focus_force()
        
            methasteprog_frame = Frame(self.methasteprog_window)
            self.methasteprog = Progressbar(methasteprog_frame,orient=HORIZONTAL, length=320, mode='determinate')
            methasteprog_frame.pack()
            self.methasteprog_window.resizable(width=False,height=False)
            starttime = time.time()
            now = starttime
            evalpersec = 1
            self.methaste_go = True
            predvar = 0
            predicttime = []
            self.methasteprog.grid(row=0,column=0,columnspan=5)
            self.prog_clock = Label(methasteprog_frame,text='')
            self.prog_clock.grid(row=1,column=1,pady=2,padx=2,sticky=W)
            prog_label = Label(methasteprog_frame,text='Run Time:')
            pred_label = Label(methasteprog_frame,text='Prediction:')
            prog_label.grid(row=1,column=0,pady=2,padx=2,sticky=E)
            pred_label.grid(row=1,column=2,pady=2,padx=2,sticky=E)
            
            
            self.pred_clock = Label(methasteprog_frame,text='')
            self.pred_clock.grid(row=1,column=3,pady=2,padx=2,sticky=W)
            self.eval_amount = Label(methasteprog_frame,text='0/'+str(Nruns))
            self.eval_amount.grid(row=1,column=4,pady=2,padx=2,sticky=W)       
            self.methasteprog["maximum"]=Nruns
            self.methasteprog.grab_set()
            self.methastecounter = 0
            
            for i in range(0,Nruns):
                pre = time.time()
                #For the range of Nruns run an SSA and record it
                X_Array = self.GeneModel.SSA(x0,t,parameters,case,TimeVar)
                S1[:,i] = X_Array[0,:]
                S2[:,i] = X_Array[1,:]
                S3[:,i] = X_Array[2,:]
                R[:,i] = X_Array[3,:]
                
                self.methasteprog.step(amount=1)  #update prog bar
                self.methasteprog_window.update()  
               # tot_param = np.zeros((1,8))
                #accept_param = np.zeros((1,8))
                self.methastecounter += 1
                
                nowold = now
                
                now = int(np.floor(time.time()-starttime))
                evaltime = time.time() - pre
               
                predtime = (Nruns)*evaltime
                predicttime.append(predtime)
                if predtime != []:
                    predav = sum(predicttime)/len(predicttime)
                
                m, s = divmod(now, 60)
                h, m = divmod(m, 60)
                curr = "%d:%02d:%02d" % (h, m, s)
                self.prog_clock.configure(text=curr)
                
                m, s = divmod(predav, 60)
                h, m = divmod(m, 60)
                pred = "%d:%02d:%02d" % (h, m, s)            
                self.pred_clock.configure(text=pred)
               
                self.eval_amount.configure(text=str(self.methastecounter)+'/'+str(Nruns))
                    

            self.methasteprog_window.destroy()        #destroy the GUI
            self.focus_force() #reurn focus to GUI

        
            f = tkFileDialog.asksaveasfile(mode='w',defaultextension='.csv')
            
            
            a = R.T
            
            
            index = np.linspace(float(self.x0_t0.get()),float(self.x0_tf.get())-1.0,ind)
            index = np.floor(index)
            
            a2 =  a[:,index.astype(int)]
  
            
            if f != None:
                np.savetxt(f,a2,delimiter=',',fmt='%1.4f',newline='\r\n')    
    
    def About(self):        
        self.Helpbox = Toplevel(self)
        self.Helpbox.geometry("710x235+300+10")
        self.Helpbox.title('')
        self.Helpbox.iconbitmap('qbio_icon_9nD_icon.ico')
        self.HelpboxFrame = Frame(self.Helpbox,relief=RAISED)
        self.HelpboxFrame.pack(fill=BOTH)
        
        self.helpcanvas = Canvas(self.HelpboxFrame,width=1200,height=900)
        self.helpcanvas.pack(anchor=NW,fill=BOTH)
        self.aboutgif = PhotoImage(file='about.gif')   
        self.helpcanvas.image =self.aboutgif
        image_on_canvas = self.helpcanvas.create_image(0, 0,anchor=NW, image=self.aboutgif)        
        
    def Help(self):
        self.Helpbox = Toplevel(self)
        self.Helpbox.geometry("1200x900+300+10")
        self.Helpbox.title('Help')
        self.Helpbox.iconbitmap('qbio_icon_9nD_icon.ico')
        self.HelpboxFrame = Frame(self.Helpbox,relief=RAISED)
        self.HelpboxFrame.pack(fill=BOTH)
        
        self.helpcanvas = Canvas(self.HelpboxFrame,width=1200,height=900)
        self.helpcanvas.pack(anchor=NW,fill=BOTH)
        self.helpgif = PhotoImage(file='test.gif')   
        self.helpcanvas.image =self.helpgif
        image_on_canvas = self.helpcanvas.create_image(0, 0,anchor=NW, image=self.helpgif)


    def Info(self):
        self.Infobox = Toplevel(self)
        self.Infobox.geometry("1270x550+300+10")
        self.Infobox.title('Model Info')
        self.Infobox.iconbitmap('qbio_icon_9nD_icon.ico')
        self.InfoboxFrame = Frame(self.Infobox,relief=RAISED)
        self.InfoboxFrame.pack(fill=BOTH)
        
        
        self.infocanvas = Canvas(self.InfoboxFrame,height=550,width=1260)
        self.infocanvas.pack(anchor=NW,fill=BOTH)
        self.infogif = PhotoImage(file='help_gui.gif')   
        self.infocanvas.image =self.infogif
        image_on_canvas = self.infocanvas.create_image(0, 0,anchor=NW, image=self.infogif)

    
    def Save_Config(self):
        f = open(os.path.normpath(os.path.dirname(os.path.abspath(__file__))+'\GUI.config.txt'),'w')
        try:
            f.write("filepath ="+self.filename+"\n")
        except:
            f.write("filepath =")
        f.write("S1 = "+self.x0_S1.get()+"\n")
        f.write("S2 = "+self.x0_S2.get()+"\n")
        f.write("S3 = "+self.x0_S3.get()+"\n")
        f.write("Rna = "+self.x0_R.get()+"\n")
        f.write("t0 = " + self.x0_t0.get() + "\n")
        f.write("tf = " + self.x0_tf.get() + "\n")    
        f.write("input = " + self.Input_entry.get() + "\n")
        f.write("k12_e = " + self.k12_entry.get() + "\n")
        f.write("k23_e = " + self.k23_entry.get() + "\n")
        f.write("k21_e = " + self.k21_entry.get() + "\n")
        f.write("k32_e = " + self.k32_entry.get() + "\n")
        f.write("kr2_e = " + self.kr2_entry.get() + "\n")
        f.write("kr3_e = " + self.kr3_entry.get() + "\n")
        f.write("gamma = " + self.gamma_entry.get() + "\n")
        f.write("beta = " + self.beta_entry.get() + "\n")
        f.write("mink12 = " + str(self.bounds[0,0]) + "\n")
        f.write("mink23 = " + str(self.bounds[0,1]) + "\n")
        f.write("mink21 = " + str(self.bounds[0,2]) + "\n")
        f.write("mink32 = " + str(self.bounds[0,3])  + "\n")
        f.write("minkr2 = " + str(self.bounds[0,4])  + "\n")
        f.write("minkr3 = " + str(self.bounds[0,5])  + "\n")
        f.write("miny = " + str(self.bounds[0,6])  + "\n")
        f.write("minb = " + str(self.bounds[0,7] ) + "\n")
        
        f.write("maxk12 = " + str(self.bounds[1,0]) + "\n")
        f.write("maxk23 = " + str(self.bounds[1,1] )+ "\n")
        f.write("maxk21 = " + str(self.bounds[1,2] )+ "\n")
        f.write("maxk32 = " + str(self.bounds[1,3] )+ "\n")
        f.write("maxkr2 = " + str(self.bounds[1,4] ) + "\n")
        f.write("maxkr3 = " + str(self.bounds[1,5])  + "\n")
        f.write("maxy = " + str(self.bounds[1,6] ) + "\n")
        f.write("maxb = " + str(self.bounds[1,7]  )+ "\n")
        f.write("helpon = " + str(self.helptips.get()) + "\n")
    

            
    ########################################################################
    # Export parameters to a proprietary data type 
    # saved as a .dat 
    ########################################################################
        
    def Export_Parameters(self):
        
        FILEOPENOPTIONS = dict(defaultextension='.dat',
                  filetypes=[('Dat file','*.dat'),('Pickle file','*.p'),('CSV','*.csv')])
                  
        filename = asksaveasfilename(parent=self.OpenFile,**FILEOPENOPTIONS)
        if filename == "":
            return
       
        k12 = self.k12_entry.get()
        k23 = self.k23_entry.get()
        k21 = self.k21_entry.get()
        k32 = self.k32_entry.get()
        kr2 = self.kr2_entry.get()
        kr3 = self.kr3_entry.get()
        gamma =  self.gamma_entry.get()
        beta = self.beta_entry.get()
        
        parameters = np.array([float(k12),   float(k21),   float(k23),   float(k32), float(kr2), float(kr3),  float(gamma),float(beta)])

        x0 = [float(self.x0_S1.get()), float(self.x0_S2.get()), float(self.x0_S3.get()),float(self.x0_R.get())]
        InputFun = self.Input_entry.get()
        t = [float(self.x0_t0.get()),float(self.x0_tf.get())]
        if filename[filename.index('.'):] == ".p"   :    
            f = open(filename,"wb")
            pickle.dump(parameters,f)
            pickle.dump(x0,f)
            pickle.dump(InputFun,f)
            pickle.dump(t,f)
            f.close
            
        if filename[filename.index('.'):] == ".csv"   :    
            f = open(filename,"wb")
            parstring = ''.join(str(e)+',' for e in parameters.tolist())
            f.write('parameters,'+parstring+'\r\n')
            x0string = ''.join(str(e)+',' for e in x0)
            f.write('x0,'+x0string+'\r\n')
            f.write('input function,' + InputFun+'\r\n')
            tstring = ''.join(str(e)+',' for e in t)
            f.write('time,'+tstring+'\r\n')
            f.close
            
                
        if filename[filename.index('.'):] == ".dat":
        
            #define the data type, numpy can do this automatically
            self.dt = np.dtype([('parameters',[('k12', float), ('k21', float), ('k23', float), ('k32', float), ('kr2', float), ('kr3', float), ('gamma', float), ('bet', float)]
                                ),
                                ('x0',[('S1',float), ('S2', float),('S3', float), ('R', float)]),
                                ('time',[('t0',float),('tf',float)]),
                                ('input',[('fun',np.dtype('|S33'))])])
            #filename, ask the user for the file via Save As dialogue from tkinter

            x = np.zeros((1,), dtype=self.dt)
            # enter the current parameters into the new data type structure
            x['parameters']['k12'] = float(self.k12_entry.get())
            x['parameters']['k21'] = float(self.k21_entry.get())
            x['parameters']['k23'] = float(self.k23_entry.get())
            x['parameters']['k32'] = float(self.k32_entry.get())
            x['parameters']['kr2'] = float(self.kr2_entry.get())
            x['parameters']['kr3'] = float(self.kr3_entry.get())
            
            x['parameters']['bet'] = float(self.beta_entry.get())
            x['parameters']['gamma'] = float(self.gamma_entry.get())
            
            x['x0']['S1'] = float(self.x0_S1.get())
            x['x0']['S2'] = float(self.x0_S2.get())
            x['x0']['S3'] = float(self.x0_S3.get())
            x['x0']['R'] = float(self.x0_R.get())
            
            x['time']['t0'] = float(self.x0_t0.get())
            x['time']['tf'] = float(self.x0_tf.get())
            
            
            x['input']['fun'] =self.Input_entry.get()
    
            x.tofile(filename,) #save te data to a file.dat
    
    def Clear_Par_Entries(self):
            self.k12_entry.delete(0,'end')
            self.k23_entry.delete(0,'end')
            self.k21_entry.delete(0,'end')
            self.k32_entry.delete(0,'end')
            self.kr2_entry.delete(0,'end')
            self.kr3_entry.delete(0,'end')
            self.gamma_entry.delete(0,'end') 
            self.beta_entry.delete(0,'end')      
 
 
    def Clear_All_Entries(self):
            self.k12_entry.delete(0,'end')
            self.k23_entry.delete(0,'end')
            self.k21_entry.delete(0,'end')
            self.k32_entry.delete(0,'end')
            self.kr2_entry.delete(0,'end')
            self.kr3_entry.delete(0,'end')
            self.gamma_entry.delete(0,'end') 
            self.beta_entry.delete(0,'end')        
            
            self.x0_S1.delete(0,'end')    
            self.x0_S2.delete(0,'end')    
            self.x0_S3.delete(0,'end')    
            self.x0_R.delete(0,'end')    

    
            self.x0_t0.delete(0,'end')
            self.x0_tf.delete(0,'end') 
            self.Input_entry.delete(0,'end') 
    ########################################################################
    # Import parameters from a proprietary data type 
    # saved as a .dat 
    ########################################################################
        
    def Import_Parameters(self):
        #define the data type (save memory if the user doesnt call this function or export)
        self.dt = np.dtype([('parameters',[('k12', float), ('k21', float), ('k23', float), ('k32', float), ('kr2', float), ('kr3', float), ('gamma', float), ('bet', float)]
                            ),
                            ('x0',[('S1',float), ('S2', float),('S3', float), ('R', float)]),
                            ('time',[('t0',float),('tf',float)]),
                            ('input',[('fun',np.dtype('|S33'))])])
        #ask the file name in open dialog
        #define a default extenstion of .dat
        FILEOPENOPTIONS = dict(defaultextension='.dat',
                  filetypes=[('Dat file','*.dat'),('Pickle file','*.p')])
        # askopen for file .dat
        filename = askopenfilename(parent=self.OpenFile,**FILEOPENOPTIONS)
        if filename == "":
            return
        
        if filename[filename.index('.'):] == ".p":  
            f = open(filename,'rb')
            try:
                parameters = pickle.load(f)
                x0 = pickle.load(f)
                InputFun = pickle.load(f)
                t = pickle.load(f)
            except:
                tkMessageBox.showerror(title='Incorrect Format',message='Something went wrong with loading this pickle file (.p), please double check the file.')
                return
            self.Clear_All_Entries()
            self.k12_entry.insert(0,str(parameters[0]))
            self.k21_entry.insert(0,str(parameters[1]))
            self.k23_entry.insert(0,str(parameters[2]))
            self.k32_entry.insert(0,str(parameters[3]))
            self.kr2_entry.insert(0,str(parameters[4]))
            self.kr3_entry.insert(0,str(parameters[5]))
            self.beta_entry.insert(0,str(parameters[7]))
            self.gamma_entry.insert(0,str(parameters[6]))  
            
            self.x0_S1.insert(0,str(x0[0]))
            self.x0_S2.insert(0,str(x0[1]))
            self.x0_S3.insert(0,str(x0[2]))
            self.x0_R.insert(0,str(x0[3]))
    
            self.x0_t0.insert(0,str(t[0]))
            self.x0_tf.insert(0,str(t[1]))
       
            self.Input_entry.insert(0,InputFun)    
            
        if filename[filename.index('.'):] == ".dat":
            x = np.fromfile(filename,dtype = self.dt)
            if not x['parameters']:
                tkMessageBox.showerror(title='Incorrect Format',message='Something is wrong with this data file (.dat), please double check the file you loaded')
                return
            if not x['x0']:
                tkMessageBox.showerror(title='Incorrect Format',message='Something is wrong with this data file (.dat), please double check the file you loaded')
                return
            if not x['input']:
                tkMessageBox.showerror(title='Incorrect Format',message='Something is wrong with this data file (.dat), please double check the file you loaded')
                return
            if not x['time']:
                tkMessageBox.showerror(title='Incorrect Format',message='Something is wrong with this data file (.dat), please double check the file you loaded')
                return
                
            #insert the floats in the entrys        
            self.Clear_All_Entries()
            
            self.k12_entry.insert(0,str(x['parameters']['k12'])[1:-1])
            self.k21_entry.insert(0,str(x['parameters']['k21'])[1:-1])
            self.k23_entry.insert(0,str(x['parameters']['k23'])[1:-1])
            self.k32_entry.insert(0,str(x['parameters']['k32'])[1:-1])
            self.kr2_entry.insert(0,str(x['parameters']['kr2'])[1:-1])
            self.kr3_entry.insert(0,str(x['parameters']['kr3'])[1:-1])
            self.beta_entry.insert(0,str(x['parameters']['bet'])[1:-1])
            self.gamma_entry.insert(0,str(x['parameters']['gamma'])[1:-1])
            
            
            self.x0_S1.insert(0,str(x['x0']['S1'])[1:-1])
            self.x0_S2.insert(0,str(x['x0']['S2'])[1:-1])
            self.x0_S3.insert(0,str(x['x0']['S3'])[1:-1])
            self.x0_R.insert(0,str(x['x0']['R'])[1:-1])
    
            self.x0_t0.insert(0,str(x['time']['t0'])[1:-1])
            self.x0_tf.insert(0,str(x['time']['tf'])[1:-1])
       
            self.Input_entry.insert(0,str(x['input']['fun'])[2:-2])


    def swapimage(self): #define a function to switch model images based on model number selection
        modelnum = self.case.get() #obtain the current case number
        self.modelgif = self.modelimages[modelnum] #redefine the current image pointer
        self.modelcanvas.itemconfigure(self.image_on_canvas, image = self.modelgif) #update the image function
        self.fig.add_subplot()
        
        
        
    def activate_prediction_boxes(self): #function to gray out SSA and FSP boxes
        if self.predict.get() == 1: #reconfigure each box 1=SSA selected
            self.SSA_entry.configure(state=NORMAL)
            self.FSP_entry.configure(state=DISABLED)

        if self.predict.get() == 2: # 2 FSP selected
            self.SSA_entry.configure(state=DISABLED)
            self.FSP_entry.configure(state=NORMAL)
            if self.FSP_entry.get() == u"ε":
                self.FSP_entry.delete(0,'end')
            
        if self.predict.get() == 0: # 0 ODE selceted
            self.SSA_entry.configure(state=DISABLED)
            self.FSP_entry.configure(state=DISABLED)
    
    def activate_fitting_boxes(self): #function to gray out SSA and FSP boxes
        if self.odeorfspvar.get() == 0: #reconfigure each box 1=SSA selected
            self.FittingFSP_entry.configure(state=DISABLED)
        
        if self.odeorfspvar.get() == 1: # 2 FSP selected
            self.FittingFSP_entry.configure(state=NORMAL)
            if self.FittingFSP_entry.get() == u"ε":
                self.FittingFSP_entry.delete(0,'end')
        
            
            
    ########################################################################
    # Update Param Callback
    # clears all entries and then replaces with search parameters
    ########################################################################
                      
                        
    def Update_Parameters(self):
        fitnum = self.fittype.get()
        odeorfsp = self.Fitting.get()
        
        self.k12_entry.delete(0,'end')
        self.k23_entry.delete(0,'end')
        self.k21_entry.delete(0,'end')
        self.k32_entry.delete(0,'end')
        self.kr2_entry.delete(0,'end')
        self.kr3_entry.delete(0,'end')
        self.gamma_entry.delete(0,'end') 
        self.beta_entry.delete(0,'end')
        
        if odeorfsp ==0:
            if fitnum ==2:
                
                self.k12_entry.insert(0,self.Best_Par_MetHaste[3][0])
                self.k23_entry.insert(0,self.Best_Par_MetHaste[3][1])
                self.k21_entry.insert(0,self.Best_Par_MetHaste[3][2])
                self.k32_entry.insert(0,self.Best_Par_MetHaste[3][3])
                self.kr2_entry.insert(0,self.Best_Par_MetHaste[3][4])
                self.kr3_entry.insert(0,self.Best_Par_MetHaste[3][5])
                self.gamma_entry.insert(0,self.Best_Par_MetHaste[3][6])
                self.beta_entry.insert(0,self.Best_Par_MetHaste[3][7])
                
        if odeorfsp ==0:
            if fitnum ==1:
    
                self.k12_entry.insert(0,self.Best_Par_SA[0][0])
                self.k23_entry.insert(0,self.Best_Par_SA[0][1])
                self.k21_entry.insert(0,self.Best_Par_SA[0][2])
                self.k32_entry.insert(0,self.Best_Par_SA[0][3])
                self.kr2_entry.insert(0,self.Best_Par_SA[0][4])
                self.kr3_entry.insert(0,self.Best_Par_SA[0][5])
                self.gamma_entry.insert(0,self.Best_Par_SA[0][6]) 
                self.beta_entry.insert(0,self.Best_Par_SA[0][7])
                
        if odeorfsp ==0:
            if fitnum ==0:
                
                self.k12_entry.insert(0,self.Best_Par_ODE[0])
                self.k23_entry.insert(0,self.Best_Par_ODE[1])
                self.k21_entry.insert(0,self.Best_Par_ODE[2])
                self.k32_entry.insert(0,self.Best_Par_ODE[3])
                self.kr2_entry.insert(0,self.Best_Par_ODE[4])
                self.kr3_entry.insert(0,self.Best_Par_ODE[5])
                self.gamma_entry.insert(0,self.Best_Par_ODE[6])  
                self.beta_entry.insert(0,self.Best_Par_ODE[7])
                


    def bringup_help(self):
        manual = Toplevel()
        self.rm_popup()
        manual.title("Help Manual")
        manual.geometry('400x700')
        manual.iconbitmap('qbio_icon_9nD_icon.ico')
        manualFrame = Frame(manual,relief=RAISED)
        manualFrame.pack(fill=BOTH)
        
        model_description_text = Button(manualFrame,text = "Model Description",bd=1,highlightthickness=1,overrelief=SUNKEN)
        model_description_text.grid(column=0,row=0,padx=2,pady=2,sticky=W)
        prediction_text = Label(manualFrame,text = "Stochasitic Simulation Algroithm")
        prediction_text.grid(column=0,row=1,padx=2,pady=2,sticky=W)        
        prediction_text_1 = Label(manualFrame,text = "Finite State Projection")
        prediction_text_1.grid(column=0,row=2,padx=2,pady=2,sticky=W)                
        
        
    def rm_popup(self):

        
        self.help_menu.destroy()
        self.menuopen = False
        
    def popup_help(self):
    
        widget = self.winfo_containing(self.winfo_pointerx(),self.winfo_pointery())
        class_type = widget.winfo_class()    
        
        itt = ""
        itt2 = ""
        itt3 = ""
        itt4 = ""
        itt5 = ""
        if class_type == "Label":
            text = widget.cget("text")
            
            if text == "Model:":
                itt = " Model type, corresponding to which rate is time dependent "
            if "Input F(t)" in text:
                itt = " Function of time fed to the selected rate parameter "
            if text == "Initial Values":
                itt  = " Starting points for all the described states, typically "
                itt2 = " S1: 1, S2: 0, S3: 0, R: N "  
                
            if "S1" in text: 
                itt = " State 1: Gene is off "
            if "S2" in text: 
                itt = " State 2: Gene is on with rate kr2 "            
            if "S3" in text: 
                itt = " State 3: Gene is on with rate kr3 " 
            if text == "R" :
                itt = " mRNA amount " 
            if "t0" in text:
                itt = " Time Zero "
            if "tf" in text: 
                itt = " Final Time "
            if "Prediction" in text:
                itt  = " What type of prediction to use: Ordinary Differential Equation (ODE), "
                itt2 = " Stochastic Simulation Algorithm (SSA), or Finite State Projection (FSP) "
                
        if class_type == "TLabel":
            text = widget.cget("text")        
            if "Parameter" in text:
                if "Search" in text:
                    itt  = " Preform different parameter searches: Built in optimization with SciPy's toolbox,  "
                    itt2 = " Simulated Annealing (SA), or Metropolis-Hastings Algrothim. If Met-Haste is chosen, "
                    itt3 = " parameter space is stored and can be viewed vs the parameters checked below. "
                
            if "Min(" in text:
                itt = " Returned minimum error found from the parameter search last ran "
            
            if "Burn" in text:
                itt = " Amount of iterations burned in the Met-Haste before starting to record outcomes "
            if "Chain" in text:
                itt = " Amount of iterations to record "
            if "Thin" in text:
                itt = " Amount of mutation attempts per each recorded chain "
            if "Mut R" in text:
                itt = " The mutation chance for each parameter, 0.2 corresponds to 20% "
            
                
                
        if class_type == "Canvas":
            itt  = " 3 State Gene Bursting Model: Each square represents a an expression of a gene: "
            itt2 = " S1 = Gene is off, S2 = Gene is on at with transcription rate kr2, S3 = Gene is "
            itt3 = " on at rate kr3. Both S2 and S3 create mRNA represented by the squiggly line. "
            itt4 = " mRNA can then decay at rate Gamma. "
        
        
            
        if itt != "":
            self.help_menu = Toplevel(self.parent)
            self.help_menu.overrideredirect(1)
            s = Style()
            s.configure("help.TFrame", background = "white")
            s.configure("help.TLabel",background = "white")
            self.help_menu.geometry('+'+str(self.winfo_pointerx())+'+'+str(self.winfo_pointery()))
            help_menu_frame = Frame(self.help_menu,relief=RIDGE,style = "help.TFrame")
            self.post_x = self.winfo_pointerx()
            self.post_y = self.winfo_pointery()
            
            
            
            infotext = Label(help_menu_frame,text=itt, style = "help.TLabel")
            infotext.pack(pady=2,padx=2,anchor=W)
            if itt2 !="":
                infotext2 = Label(help_menu_frame,text=itt2,style = "help.TLabel")
                infotext2.pack(pady=2,padx=2,anchor=W)
            if itt3 !="":
                infotext3 = Label(help_menu_frame,text=itt3,style = "help.TLabel")
                infotext3.pack(pady=2,padx=2,anchor=W)
            if itt4 !="":
                infotext4 = Label(help_menu_frame,text=itt4,style = "help.TLabel")
                infotext4.pack(pady=2,padx=2,anchor=W)
            if itt5 !="":
                infotext5 = Label(help_menu_frame,text=itt5,style = "help.TLabel")
                infotext5.pack(pady=2,padx=2,anchor=W)            
                
    
            infobutton = Button(help_menu_frame,text=" More ",command=self.bringup_help,bg="White")
            
            infobutton.pack(anchor=E,pady=2,padx=2)
    
            help_menu_frame.pack(pady=1,padx=1)
            self.menuopen = True

        
    def stop_poll(self,event):
        self.offscreen = True
    def start_poll(self,event):
        self.offscreen = False
        
    def poll(self):
        
        if self.offscreen == False and self.helptips.get() == 1:
            
            if self.winfo_pointerx() == self.x and self.winfo_pointery() == self.y:
                self.count +=1
            else:
        
                if self.winfo_pointerx() != self.x:
                    self.x = self.winfo_pointerx()
                    self.count = 0
                if self.winfo_pointery() != self.y:
                    self.y = self.winfo_pointery()
                    self.count = 0
            
            if self.count >= 120 and self.menuopen == False:
                self.popup_help()
                
                
            if self.menuopen == True:
                
                if np.sqrt((self.post_x-self.winfo_pointerx())**2+(self.post_y-self.winfo_pointery())**2) >50:
                    self.rm_popup()
        self.afer_id = self.after(1,self.poll)        

                
                
    def MH_options(self):
        if self.MH_options_open == True:
            self.MH_window.focus()
            return
            
        
        
        self.MH_window = Toplevel()
        self.MH_window.protocol("WM_DELETE_WINDOW", self.MH_close)
        self.MH_window.iconbitmap('qbio_icon_9nD_icon.ico') #set icon
        self.MH_window.title("MH")
        MH_frame = Frame(self.MH_window,relief=RAISED)
        MH_frame.pack()
        
        self.MH_window.focus()        
        
        self.MH_textbox = Text(MH_frame,width=30,height=10,state=DISABLED)
        self.MH_textbox.grid(column=0,row=0,padx=2,pady=2,columnspan=4)
        
        Plot_label = Label(MH_frame,text="Parameters")
        Plot_label.grid(column=0,row=3,padx=6,pady=2,columnspan=2,sticky=W)
        
        Acc_Plot_label = Label(MH_frame,text="Autocorrelations")
        Acc_Plot_label.grid(column=2,row=3,padx=2,pady=2,columnspan=2,sticky=W)
                
        Plot_par = Menubutton(MH_frame,relief=RAISED,text="Params to Plot")
        Plot_par_menu=Menu(Plot_par,tearoff=0)
        Plot_par["menu"] = Plot_par_menu
        
        self.Plot_par_boxes =[]
        for i in range(0,8):
            self.Plot_par_boxes.append(IntVar(value=1))
        
        labels = ['k12','k23','k21','k32','kr2','kr3','gmm','beta']
        for i in range(len(self.Plot_par_boxes)):
            Plot_par_menu.add_checkbutton(label = labels[i],onvalue=1,offvalue=0,variable=self.Plot_par_boxes[i])
            
        Plot_par.grid(column=0,row=4,padx=6,pady=2,columnspan=2,sticky=W)
        
        colormode = [('Solid Color','Solid'),('Tricolor','Tri'),('Heatmap','Jet'),('Greyscale','Greys')]
        self.pcolor = StringVar()
        self.pcolor.set('Solid')       
        Plot_par_c = Menubutton(MH_frame,relief=RAISED,text="Color Options")
        Plot_par_c_menu=Menu(Plot_par_c,tearoff=0)
        Plot_par_c["menu"] = Plot_par_c_menu                   
        for text, colormodes in colormode:
            Plot_par_c_menu.add_radiobutton(label = text,variable=self.pcolor,value=colormodes)
        Plot_par_c.grid(column=0,row=5,padx=6,pady=2,columnspan=2,sticky=W)
        
        
        transp = [('.01',.01),('.1',.1),('.3',.3),('.5',.5),('.7',.7),('1',1.0)]
        self.transp_var = DoubleVar()
        self.transp_var.set(.1)       
        Plot_par_t = Menubutton(MH_frame,relief=RAISED,text="Transparency")
        Plot_par_t_menu=Menu(Plot_par_t,tearoff=0)
        Plot_par_t["menu"] = Plot_par_t_menu                   
        for text, transps in transp:
            Plot_par_t_menu.add_radiobutton(label = text,variable=self.transp_var,value=transps)
        Plot_par_t.grid(column=0,row=6,padx=6,pady=2,columnspan=2,sticky=W)       
        
        Plot_Par_space_button = Button(MH_frame,text="Plot Par Space",command = self.MH_par_space_call)
        Plot_Par_space_button.grid(column=0,row=8,padx=6,pady=6,columnspan=2,sticky=W)
        
        Plot_start_label = Label(MH_frame,text="Plot: ")
        self.Plot_start_entry = Entry(MH_frame,width=11)
        Plot_mid_label = Label(MH_frame,text="to")
        self.Plot_end_entry = Entry(MH_frame,width=11)
        
        Plot_start_label.grid(column=0,row=7,padx=2,pady=6)
        self.Plot_start_entry.grid(column=1,row=7,padx=2,pady=6,sticky=E)
        Plot_mid_label.grid(column=2,row=7,padx=2,pady=6)
        self.Plot_end_entry.grid(column=3,row=7,padx=2,pady=6,sticky=W)
        
        Plot_feval = Button(MH_frame,text="Plot feval",command=self.plot_MH_feval)
        Plot_feval.grid(column=2,row=8,padx=2,columnspan=2,sticky=W)
        
        Plot_Param_autocorr = Button(MH_frame,text="Parameter Accs",command=self.plot_MH_par_acc)
        Plot_Param_autocorr.grid(column=2,row=4,padx=2,columnspan=2,sticky=W)
        
        Plot_Feval_autocorr = Button(MH_frame,text="Feval Accs",command=self.plot_MH_feval_acc)
        Plot_Feval_autocorr.grid(column=2,padx=2,row=5,columnspan=2,sticky=W)
        
        
        
        Load_MH = Button(MH_frame,text="Load MH run",command=self.Load_MH)
        Load_MH.grid(column=0,row=1,padx=6,pady=2,columnspan=2,sticky=W)
        Save_MH = Button(MH_frame,text="Save MH run",command = self.Save_MH)
        Save_MH.grid(column=2,row=1,padx=6,pady=2,columnspan=2,sticky=W)
        
        Send_par = Button(MH_frame,text="Send Par to Main Window",command = self.Send_mh_par_to_main)
        Send_par.grid(column=0,row=2,padx=6,pady=2,columnspan=4,sticky=W)
        
        self.MH_options_open = True
        try: 
            self.Loaded_MH = self.Best_Par_MetHaste
            self.update_mh_options()
        except:
            self.MH_textbox.config(state=NORMAL)
            self.MH_textbox.insert(END,' No met-haste run loaded\n')
            self.MH_textbox.insert(END,' Load or Run a MH search\n')
            self.MH_textbox.config(state=DISABLED)
           
            pass
        
            
    def Save_MH(self):
        try:
            filename = tkFileDialog.asksaveasfilename(defaultextension='.csv')
            tosave = np.vstack((self.Loaded_MH[1],self.Loaded_MH[0]))
            np.savetxt(filename,tosave, fmt='%.5f',delimiter=',',newline='\n')
        except:
            return
        if filename =='':
            return            
            
    def MH_close(self):
        self.MH_options_open = False
        self.MH_window.destroy()
        
        
        
    def Load_MH(self):
        try:
            filename = tkFileDialog.askopenfilename()
        except:
            return
        
        if filename =='':
            return
        
        Pars = np.loadtxt(filename,delimiter=',')

        self.Loaded_MH = Pars[8], Pars[0:8], min(Pars[8]), Pars[0:8].T[Pars[8].tolist().index(min(Pars[8]))]
        self.update_mh_options()
        self.MH_window.focus()       
    
    def Send_mh_par_to_main(self):
        try:
            
            if self.Loaded_MH =="":
                return
            
            self.k12_entry.delete(0,'end')
            self.k23_entry.delete(0,'end')
            self.k21_entry.delete(0,'end')
            self.k32_entry.delete(0,'end')
            self.kr2_entry.delete(0,'end')
            self.kr3_entry.delete(0,'end')
            self.gamma_entry.delete(0,'end') 
            self.beta_entry.delete(0,'end')
            
            self.k12_entry.insert(0,self.Loaded_MH[3][0])
            self.k23_entry.insert(0,self.Loaded_MH[3][1])
            self.k21_entry.insert(0,self.Loaded_MH[3][2])
            self.k32_entry.insert(0,self.Loaded_MH[3][3])
            self.kr2_entry.insert(0,self.Loaded_MH[3][4])
            self.kr3_entry.insert(0,self.Loaded_MH[3][5])
            self.gamma_entry.insert(0,self.Loaded_MH[3][6])
            self.beta_entry.insert(0,self.Loaded_MH[3][7])
        except:
            return
                            


    def update_mh_options(self):
        self.MH_textbox.config(state=NORMAL)
        self.MH_textbox.delete(0.0,END)
        self.MH_textbox.insert(END,' Best Parameters\n')
        
        self.MH_textbox.insert(END,' k12 = ')
        self.MH_textbox.insert(END,self.Loaded_MH[3][0])
        self.MH_textbox.insert(END,' \n')

        self.MH_textbox.insert(END,' k21 = ')
        self.MH_textbox.insert(END,self.Loaded_MH[3][1])
        self.MH_textbox.insert(END,' \n')

        self.MH_textbox.insert(END,' k23 = ')
        self.MH_textbox.insert(END,self.Loaded_MH[3][2])
        self.MH_textbox.insert(END,' \n')

        self.MH_textbox.insert(END,' k32 = ')
        self.MH_textbox.insert(END,self.Loaded_MH[3][3])
        self.MH_textbox.insert(END,' \n')

        self.MH_textbox.insert(END,' kr2 = ')
        self.MH_textbox.insert(END,self.Loaded_MH[3][4])
        self.MH_textbox.insert(END,'\n')            

        self.MH_textbox.insert(END,' kr3 = ')
        self.MH_textbox.insert(END,self.Loaded_MH[3][5])
        self.MH_textbox.insert(END,'\n')

        self.MH_textbox.insert(END,' gamma = ')
        self.MH_textbox.insert(END,self.Loaded_MH[3][6])
        self.MH_textbox.insert(END,'\n')
        
        self.MH_textbox.insert(END,' beta = ')
        self.MH_textbox.insert(END,self.Loaded_MH[3][7])
        self.MH_textbox.insert(END,'\n')
        
        
        self.MH_textbox.insert(END,' Feval =  ')
        self.MH_textbox.insert(END,self.Loaded_MH[2])
        self.MH_textbox.insert(END,'\n')
        self.MH_textbox.config(state=DISABLED)

        self.Plot_start_entry.delete(0,'end')
        self.Plot_start_entry.insert(0,str(0))
        self.Plot_end_entry.delete(0,'end')
        self.Plot_end_entry.insert(0,str(len(self.Loaded_MH[0])))         
    
        
        

    def plot_MH_feval_acc(self):
        try:
            self.Loaded_MH[0]
        except:
            return
            
        if self.test_plot_range_mh() == False:
            tkMessageBox.showerror(title='Input Error',message='Please input a valid number for the Met-Haste plotting range.')
            return
        else:
            in1,in2 = self.test_plot_range_mh()
            
        self.fig.clf()
        self.subplot = self.fig.add_subplot(111)
        
        a = self.get_acc2(self.Loaded_MH[0][in1:in2])-np.average(self.Loaded_MH[0][in1:in2])
        self.subplot.plot(np.linspace(in1,in2,in2-in1),a,label='Feval')
        self.subplot.plot([in1,in2],[0,0],color='r')
        self.subplot.set_xlabel('Feval Autocorrelation')
        self.figurecanvas.show()
        self.fig.canvas.draw()
        self.currplot.set(2)            
        
    def test_plot_range_mh(self):
        try:
            in1 = int(self.Plot_start_entry.get())
            in2 = int(self.Plot_end_entry.get())
            if in2 > len(self.Loaded_MH[0]):
                in2 = len(self.Loaded_MH[0])
                self.Plot_end_entry.delete(0,'end')
                self.Plot_end_entry.insert(0,str(in2))
                
            if in1 <0:
                in1=0
                self.Plot_start_entry.delete(0,'end')
                self.Plot_start_entry.insert(0,str(in1))
            return in1,in2
        except:
            return False
  
 
    def plot_MH_par_acc(self):
        parchain = self.Loaded_MH[1]
      
        self.fig.clf()
        labels = ['k12','k23','k21','k32','kr2','kr3','gmm','beta']
        if self.test_plot_range_mh() == False:
            tkMessageBox.showerror(title='Input Error',message='Please input a valid number for the Met-Haste plotting range.')
            return
        else:
            in1,in2 = self.test_plot_range_mh()
            
        for i in range(0,8):
            num = '42'+str(i+1)
            self.subplot = self.fig.add_subplot(int(num))
            
            a = self.get_acc2(parchain[i][in1:in2])-np.average(parchain[i][in1:in2])
            
            
            
            
            self.subplot.plot(np.linspace(in1,in2,in2-in1)  , a)
            self.subplot.plot([in1,in2],[0,0],color='r')
            self.subplot.set_title(labels[i])
        self.figurecanvas.show()
        self.fig.tight_layout()
        self.fig.canvas.draw()
        self.currplot.set(2)                 
            

    def plot_MH_feval(self):
        
        try:
            self.Loaded_MH[0]
        except:
            return
            
        if self.test_plot_range_mh() == False:
            tkMessageBox.showerror(title='Input Error',message='Please input a valid number for the Met-Haste plotting range.')
            return
        else:
            in1,in2 = self.test_plot_range_mh()
            
        self.fig.clf()
        self.subplot = self.fig.add_subplot(111)
        
        self.subplot.plot(np.linspace(in1,in2,in2-in1)  ,self.Loaded_MH[0][in1:in2],label='Feval')
        
        self.subplot.set_xlabel('Function Iterations')
        self.figurecanvas.show()
        self.fig.canvas.draw()
        self.fig.tight_layout()
        self.currplot.set(2)            



    def MH_par_space_call(self):

        if self.test_plot_range_mh() == False:
            tkMessageBox.showerror(title='Input Error',message='Please input a valid number for the Met-Haste plotting range.')
            return
        else:
            in1,in2 = self.test_plot_range_mh()
        
        parchain = self.Loaded_MH[1]
        
        
        pco = self.pcolor.get()
        al = self.transp_var.get()      
        
        
        k12plot = self.Plot_par_boxes[0].get()
        k21plot = self.Plot_par_boxes[1].get()
        k23plot = self.Plot_par_boxes[2].get()
        k32plot = self.Plot_par_boxes[3].get()
        kr2plot = self.Plot_par_boxes[4].get()
        kr3plot = self.Plot_par_boxes[5].get()
        gmmplot = self.Plot_par_boxes[6].get()
        betplot = self.Plot_par_boxes[7].get()
        
        k12list = []
        k21list = []
        k23list = []
        k32list = []
        kr2list = []
        kr3list = []
        gmmlist = []
        betlist = []
        
        
        for i in range(in1,in2):
            if k12plot ==1:
                k12list.append(parchain[0][i])
            if k23plot ==1:
                k23list.append(parchain[1][i])
            if k21plot ==1:
                k21list.append(parchain[2][i])
            if k32plot ==1:
                k32list.append(parchain[3][i])
            if kr2plot ==1:
                kr2list.append(parchain[4][i])
            if kr3plot ==1:
                kr3list.append(parchain[5][i])
            if gmmplot ==1:
                gmmlist.append(parchain[6][i])
            if betplot ==1:
                betlist.append(parchain[7][i])
            
        
        self.fig.clf() 
        self.currplot.set(2)
        self.ParSpace_Plot(self,k12plot,k21plot,k23plot,k32plot,kr2plot,kr3plot,gmmplot,betplot,k12list,k23list,k21list,k32list,kr2list,kr3list,gmmlist,betlist,alpha =al, color =pco )             
                            
                    
                                                        
        
        
        
            
    def get_acc2(self,data,trunc=False):
        N = len(data)
        fvi = np.fft.fft(data,n=2*N)
        acf = fvi*np.conjugate(fvi)
        acf = np.fft.ifft(acf)
        acf = np.real(acf[:N])/float(N)
        if trunc:
            acf[acf<0]=0
            for i in range(1,len(acf)):
                if acf[i] > acf[i-1]:
                    acf[i] = acf[i-1]
        return acf
                
                
                
                
                
                
                
                
    def PlotInput(self):
        
        checkfun = self.test_input_function(self.Input_entry.get())
        TimeVar = self.GeneModel.TimeVar_Create(self.Input_entry.get())
        if checkfun == 0:
            if self.holdon.get() == 0:
                self.fig.clf()
            self.subplot = self.fig.add_subplot(111)
            t = np.linspace(float(self.x0_t0.get()),float(self.x0_tf.get()),10*(float(self.x0_tf.get())-float(self.x0_t0.get())))
            Timeplot=[]
            for i in range(0,len(t)):
                Timeplot.append(TimeVar(t[i]))           
            self.subplot.plot(t,Timeplot)
            
            self.subplot.set_xlabel('$Time$')
            self.fig.canvas.draw()
            self.figurecanvas.show()
            
            
    ########################################################################
    # Main GUI frame and elements
    ########################################################################

    def initUI(self):
        #define the title
        self.parent.title("qbio GUI")
        self.pack(fill=BOTH, expand=True)
        
        #initate the model variable
        self.case = IntVar()
        self.case.set(0)
        
        self.currplot = IntVar()
        self.currplot.set(0)
        
        #inticate images for the model display
        #will be encoded strings when finalized?
        self.modelimages =   []
        self.modelimages.append(PhotoImage(file="updatedmodelZero.gif"))
        self.modelimages.append(PhotoImage(file="updatedmodelOne.gif"))
        self.modelimages.append(PhotoImage(file="updatedmodelTwo.gif"))
        self.modelimages.append(PhotoImage(file="updatedmodelThree.gif"))
        self.modelimages.append(PhotoImage(file="updatedmodelFour.gif"))
#_____________________________________________________________________________        
#        #set up the model display frame
#        
        #define the top menu bar for file, fitting, prediction, figure, and help
        self.menubar = Menu(self)
        menu = Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="File", menu=menu)
        
        menu.add_command(label="Import Parameters",command=self.Import_Parameters)
        menu.add_command(label="Export Parameters",command=self.Export_Parameters)
        menu.add_command(label="Save Configuration", command=self.Save_Config)
        menu.add_separator()
        menu.add_command(label="Exit",command=self.quit)

        menu = Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="Prediction", menu=menu)
        menu.add_command(label="Simulate Data via SSA",command = self.Sim_SSA)
        menu.add_command(label="Clear Parameters",command = self.Clear_Par_Entries)
        self.master.config(menu=self.menubar)
        
        menu = Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="Fitting", menu=menu)
        menu.add_command(label="Change Bounds",command=self.Change_Bounds)
        menu.add_command(label="More MH options",command=self.MH_options)
        self.master.config(menu=self.menubar)
        
        menu = Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="Figure", menu=menu)
        menu.add_command(label="Change Plot Colormap",command=self.Change_cmap)
        
        menu.add_command(label="Export Plot to txt",command=self.Export_PlotVar)
        self.master.config(menu=self.menubar)
        
        menu = Menu(self.menubar, tearoff=0)
        self.menubar.add_cascade(label="Help", menu=menu)
        self.helptips = IntVar()
        self.helptips.set(1)
        menu.add_checkbutton(label="Help Tips",onvalue = 1, offvalue = 0,variable=self.helptips)
        
        
        menu.add_command(label="Model Info",command=self.Info)
        menu.add_command(label="About",command=self.About)
        self.master.config(menu=self.menubar)
        
        #define the bottom frame and open figure buttons
        LoadingandFigFrame = Frame(self,relief=RAISED)
        LoadingandFigFrame.grid(row=11,column=0,sticky=E,columnspan=3)
        OpenFigure = Button(LoadingandFigFrame,text="Open Figure",command=self.Open_Fig)
        OpenFigure.grid(row=0,column=5,padx=4,pady=3,columnspan=1)
        
        PlotInput = Button(LoadingandFigFrame, text = "Plot Input Function",command = self.PlotInput)
        PlotInput.grid(row=0,column=1,padx=4,pady=3)
        spacelabel = Label(LoadingandFigFrame,text=' ',width=19)
        spacelabel.grid(column=2,row=0)
        self.holdoncheck = Checkbutton(LoadingandFigFrame,text="Hold on? ",variable=self.holdon)
        self.holdoncheck.grid(row=0,column=4,padx=5,pady=2)
        self.ModelPictureFrame = Frame(self,relief = RAISED)
        self.ModelPictureFrame.grid(row=1,column=0,columnspan=2,sticky=E)
        
        #file controlling which picture is displayed for the model
            
        self.modelpictureframes(self.ModelPictureFrame)

        #import all the frames and buttons associated, stored in seperate files

        ModelsFrame= Frame(self,relief=RAISED)
        ModelsFrame.grid(row=0,column=0,sticky=W)
        self.InputFrame = Frame(self,relief=RAISED)
        self.InputFrame.grid(row=0,column=1,sticky=E)      
        self.ModelsInputsFrames(ModelsFrame,self.InputFrame)
        
        self.variablesFrame = Frame(self,relief=RAISED)
        self.variablesFrame.grid(row=3,column=1,sticky=E)
        self.VariablesFrames(self.variablesFrame)
 
        self.x0Frame = Frame(self,relief=RAISED)
        self.x0Frame.grid(row=2,column=0,columnspan=2,sticky=W)       
        self.x0Frames(self.x0Frame)
        
        PredictionsFrame= Frame(self,relief=RAISED)
        PredictionsFrame.grid(row=3,column=0,sticky=W,pady=1)
        self.PredictionFrames(PredictionsFrame)
                
        self.odeorfspvar = IntVar()
        self.odeorfspvar.set(0)
        FittingFrame=Frame(self,relief=RAISED)
        FittingFrame.grid(row=4,column=0,columnspan=2,sticky=W)        
        FittingFrame2=Frame(self,relief=RAISED)
        FittingFrame2.grid(column=1,row=4,columnspan=1,sticky=E)
        self.FittingsFrames(FittingFrame,FittingFrame2)
        plt.plot
        
#        import typeoffitframe #ode or fsp radiobutton
#        TypeofFitFrame = Frame(self,relief=RAISED)
#        TypeofFitFrame.place(x=1,y=508)
#        typeoffitframe.TypeOfFittingFrames(self,TypeofFitFrame)
#        

        self.fittype = IntVar()
        self.fittype.set(0)        
        FminFrame = Frame(self)
        FminFrame.grid(row=5,column=0,columnspan=7,pady=1,sticky=W)
        TypeFitLabel = Label(FminFrame,text="Parameter Search")
        TypeFitLabel.grid(row=0,column=0,padx=2,pady=2,sticky=W)
        s = Style()
        s.configure('My.TFrame',background = "#DCDCDC")
        s.configure('My.TFrame_2',background = "white")
        
        s.configure('My.TLabel',background = "#DCDCDC")
        s.configure('My.TEntry',bg="#DCDCDC")
        FminFrame2 = Frame(self,relief=SUNKEN,style = 'My.TFrame')
        FminFrame2.grid(row=6,column=0,columnspan=3,pady=1,sticky=E)
        OptimizationFmin = Radiobutton(FminFrame2,text="scipy.minimize    ",variable = self.fittype, value =0,command=self.activate_fittype_box,bg="#DCDCDC")
        OptimizationFmin.grid(row=1,column=0,padx=2,pady=2,columnspan = 2)
        self.fminmethod = StringVar()
        self.Fmin_Method = OptionMenu(FminFrame2,self.fminmethod,"L-BFGS-B", "TNC","SLSQP") 
        self.Fmin_Method["menu"].config(bg="#DCDCDC")
        self.Fmin_Method.config(width=8)
        self.Fmin_Method["bg"] = "#DCDCDC"
        
        self.fminmethod.set('L-BFGS-B')
        self.Fmin_Method.grid(row=1,column=2,padx=3,pady=2,columnspan = 3)
        self.MinError = Text(FminFrame2,width = 17,height = 1,bg="#F0F0F0")
        MinErrorLabel = Label(FminFrame2,text="Min(ε)",style="My.TLabel")
        MinErrorLabel.grid(row=1,column=7,padx=3,pady=2,sticky=E)
        self.MinError.grid(row=1,column=8,padx=2,pady=2,columnspan=2,sticky=W)
        self.MinError.config(state=DISABLED)

        
        
        SAFrame = Frame(self,relief=SUNKEN,style='My.TFrame')
        SAFrame.grid(row=7,column=0,columnspan=3,sticky=E)

        
        
        OptimizationSA = Radiobutton(SAFrame,text="SA              ",variable = self.fittype, value =1,command=self.activate_fittype_box, bg="#DCDCDC")
        OptimizationSA.grid(row=0,column=0,padx=2,pady=2)
        SA_label = Label(SAFrame,text="Max Iterations",style="My.TLabel")
        SA_label.grid(row=0,column=1,columnspan=3)
        self.SA_limit = Entry(SAFrame,width=7,justify=RIGHT,style="My.TEntry")
        self.SA_limit.grid(row=0,column=4,padx=3,pady=2,columnspan = 2,sticky=W)
        self.SA_limit.insert(0,"500")
        self.SA_limit.config(state=DISABLED)
        MHFrame = Frame(self,relief=SUNKEN,style='My.TFrame')
        MHFrame.grid(column=0,row=8,columnspan=3,sticky=E)
        
        OptimizationMetHaste = Radiobutton(MHFrame,text="Met-Haste",variable = self.fittype, value =2,command=self.activate_fittype_box,bg="#DCDCDC")
        OptimizationMetHaste.grid(row=3,column=0,padx=2,pady=2)
        
        self.burn_label = Label(MHFrame, text="Burn",style="My.TLabel")
        self.burn_label.grid(row=3,column=1,padx=2,pady=2)
        
        self.chain_label = Label(MHFrame, text="Chain",style="My.TLabel")
        self.chain_label.grid(row=3,column=3,padx=1,pady=2)
        
        self.thin_label = Label(MHFrame, text="Thin",style="My.TLabel")
        self.thin_label.grid(row=3,column=5,padx=2,pady=2)
        
        self.thin_label = Label(MHFrame, text="Mut R",style="My.TLabel")
        self.thin_label.grid(row=3,column=7,padx=2,pady=2)        
        
        self.N_Burn = Entry(MHFrame,width=5,justify=RIGHT,state=NORMAL)
        self.N_Burn.grid(row=3,column=2,padx=4,pady=2)
        
        self.N_Chain = Entry(MHFrame,width=5,justify=RIGHT,state=NORMAL)
        self.N_Chain.grid(row=3,column=4,padx=4,pady=2)
        
        self.N_Thin = Entry(MHFrame,width=5,justify=RIGHT,state=NORMAL)
        self.N_Thin.grid(row=3,column=6,padx=4,pady=2)
        
        self.Mut_Rate = Entry(MHFrame,width=5,justify=RIGHT,state=NORMAL)
        self.Mut_Rate.grid(row=3,column=8,padx=6,pady=2)        
        
        SpaceLabel = Label(SAFrame,text='     ',width = 30,style="My.TLabel")
        SpaceLabel.grid(row=0,column=9,padx=3,pady=2)       
        
        self.N_Burn.insert(0,"10")   
        self.N_Thin.insert(0,"10")
        self.N_Chain.insert(0,"100") 
        self.Mut_Rate.insert(0,".2")
        
        self.N_Burn.config(state=DISABLED)   
        self.N_Thin.config(state=DISABLED)
        self.N_Chain.config(state=DISABLED)   
        self.Mut_Rate.config(state=DISABLED)  
        
        
        self.MetHasteParSpace = Button(MHFrame,text="View Parameter Space",state=DISABLED,command=self.PlotParameterSpace,bg="#DCDCDC")
        self.MetHasteParSpace.grid(row=4,column=5,padx=3,pady=3,columnspan=4)
        
        checkboxframe = Frame(MHFrame,style="My.TFrame")
        checkboxframe.grid(row=4,column=0,columnspan=5,padx=2,pady=2)
        self.k12checkvar = IntVar()
        self.k12checkvar.set(0)
        self.k12check = Checkbutton(checkboxframe,text="k12",variable=self.k12checkvar,bg="#DCDCDC")
        self.k12check.grid(row=0,column=0,padx=2,pady=2)

        self.k23checkvar = IntVar()
        self.k23checkvar.set(0)
        self.k23check = Checkbutton(checkboxframe,text="k23",variable=self.k23checkvar,bg="#DCDCDC")
        self.k23check.grid(row=0,column=1,padx=2,pady=2)
        
        self.k21checkvar = IntVar()
        self.k21checkvar.set(0)
        self.k21check = Checkbutton(checkboxframe,text="k21",variable=self.k21checkvar,bg="#DCDCDC")
        self.k21check.grid(row=0,column=2,padx=2,pady=2)

        self.k32checkvar = IntVar()
        self.k32checkvar.set(0)
        self.k32check = Checkbutton(checkboxframe,text="k32",variable=self.k32checkvar,bg="#DCDCDC")
        self.k32check.grid(row=1,column=0,padx=2,pady=2)
    
        self.kr2checkvar = IntVar()
        self.kr2checkvar.set(0)
        self.kr2check = Checkbutton(checkboxframe,text="kr2",variable=self.kr2checkvar,bg="#DCDCDC")
        self.kr2check.grid(row=1,column=1,padx=2,pady=2)

        self.kr3checkvar = IntVar()
        self.kr3checkvar.set(0)
        self.kr3check = Checkbutton(checkboxframe,text="kr3",variable=self.kr3checkvar,bg="#DCDCDC")
        self.kr3check.grid(row=1,column=2,padx=2,pady=2)
        
        self.gmmcheckvar = IntVar()
        self.gmmcheckvar.set(0)
        self.gmmcheck = Checkbutton(checkboxframe,text=" γ ",variable=self.gmmcheckvar,bg="#DCDCDC")
        self.gmmcheck.grid(row=1,column=3,padx=2,pady=2)

        self.betcheckvar = IntVar()
        self.betcheckvar.set(0)
        self.betcheck = Checkbutton(checkboxframe,text=" β ",variable=self.betcheckvar,bg="#DCDCDC")
        self.betcheck.grid(row=0,column=3,padx=2,pady=2)
        
        self.k12check.configure(state=DISABLED)
        self.k21check.configure(state=DISABLED)
        self.k23check.configure(state=DISABLED)
        self.k32check.configure(state=DISABLED)
        self.kr2check.configure(state=DISABLED)
        self.kr3check.configure(state=DISABLED)
        self.gmmcheck.configure(state=DISABLED)
        self.betcheck.configure(state=DISABLED)

        ButtonsFrame = Frame(self,relief=RAISED)
        ButtonsFrame.grid(row=9,column=0,columnspan=3,sticky=E)
        
        FitRun = Button(ButtonsFrame,text=" Run Search ",command=self.Run_Fit)
        FitRun.grid(row=0,column=0,padx=4,pady=3,columnspan=1,sticky=W)
        
        UpdatePar = Button(ButtonsFrame,text="Update Parameters",command=self.Update_Parameters)
        UpdatePar.grid(row=0,column=1,padx=6,pady=3,columnspan=3)
        
        SavePar = Button(ButtonsFrame,text="Save Parameters",command=self.Export_Parameters)
        SavePar.grid(row=0,column=4,padx=6,pady=3,columnspan=3)
        
        Moments = Button(ButtonsFrame,text=" 2nd Moments ",command=self.Second_Moment_Call)
        Moments.grid(row=0,column=7,padx=3,pady=3,columnspan=2)
        
        try:
            
            with open(os.path.normpath(os.path.dirname(os.path.abspath(__file__))+'\GUI.config.txt')) as config:
               self.Clear_All_Entries()
               for line in config:
                    
                    if "filepath" in line:
                        ind = line.find('=')
                        
                        filename = line[ind+1:]
                        self.filename = filename[:-1]
                        try:
                            self.mRNA_Hist = self.GeneModel.Import_mRNA(self.filename)
                            self.FileEntry.config(state=NORMAL)
                            self.FileEntry.delete(0.0,'end')
                            self.FileEntry.insert(INSERT,filename+" loaded")
                            self.FileEntry.config(state=DISABLED)
                        except:
                            x=1
                  
                    if "S1" in line:
                        ind = line.find('=')                  
                        a = line[ind+2:]
                        self.x0_S1.insert(0,float(a))
                    if "S2" in line:
                        ind = line.find('=')                  
                        a = line[ind+2:]
                        self.x0_S2.insert(0,float(a))
                    if "S3" in line:
                        ind = line.find('=')                  
                        a = line[ind+2:]
                        self.x0_S3.insert(0,float(a))
                    if "Rna" in line:
                        ind = line.find('=')                  
                        a = line[ind+2:]
                        self.x0_R.insert(0,float(a)) 
                    if "t0" in line:
                        ind = line.find('=')                  
                        a = line[ind+2:]
                        self.x0_t0.insert(0,float(a))
                    if "tf" in line:
                        ind = line.find('=')                  
                        a = line[ind+2:]
                        self.x0_tf.insert(0,float(a))
                    if "input" in line:
                        ind = line.find('=')                  
                        a = line[ind+2:]
                        self.Input_entry.insert(0,a)
                    if "k12_e" in line:
                        ind = line.find('=')                  
                        a = line[ind+2:]
                        self.k12_entry.insert(0,float(a))
                    if "k23_e" in line:
                        ind = line.find('=')                  
                        a = line[ind+2:]
                        self.k23_entry.insert(0,float(a))
                    if "k21_e" in line:
                        ind = line.find('=')                  
                        a = line[ind+2:]
                        self.k21_entry.insert(0,float(a))
                    if "kr2_e" in line:
                        ind = line.find('=')                  
                        a = line[ind+2:]
                        self.kr2_entry.insert(0,float(a))
                    if "kr3_e" in line:
                        ind = line.find('=')                  
                        a = line[ind+2:]
                        self.kr3_entry.insert(0,float(a))
                    if "k32_e" in line:
                        ind = line.find('=')                  
                        a = line[ind+2:]
                        self.k32_entry.insert(0,float(a))
                    if "gamma" in line:
                        ind = line.find('=')                  
                        a = line[ind+2:]
                        self.gamma_entry.insert(0,float(a))
                    if "beta" in line:
                        ind = line.find('=')                  
                        a = line[ind+2:]
                        self.beta_entry.insert(0,float(a))
                    if "mink12" in line:
                        ind = line.find('=')                  
                        a = line[ind+2:]
                        self.bounds[0,0] = float(a)
                    if "mink23" in line:
                        ind = line.find('=')                  
                        a = line[ind+2:]
                        self.bounds[0,1] = float(a)
                    if "mink21" in line:
                        ind = line.find('=')                  
                        a = line[ind+2:]
                        self.bounds[0,2] = float(a)
                    if "mink32" in line:
                        ind = line.find('=')                  
                        a = line[ind+2:]
                        self.bounds[0,3] = float(a)                  
                    if "minkr2" in line:
                        ind = line.find('=')                  
                        a = line[ind+2:]
                        self.bounds[0,4] = float(a)  
                    if "minkr3" in line:
                        ind = line.find('=')                  
                        a = line[ind+2:]
                        self.bounds[0,5] = float(a)    
                    if "miny" in line:
                        ind = line.find('=')                  
                        a = line[ind+2:]
                        self.bounds[0,6] = float(a)    
                    if "minb" in line:
                        ind = line.find('=')                  
                        a = line[ind+2:]
                        self.bounds[0,7] = float(a)     
                    if "maxk12" in line:
                        ind = line.find('=')                  
                        a = line[ind+2:]
                        self.bounds[1,0] = float(a)                                                   
                    if "maxk23" in line:
                        ind = line.find('=')                  
                        a = line[ind+2:]
                        self.bounds[1,1] = float(a)   
                    if "maxk21" in line:
                        ind = line.find('=')                  
                        a = line[ind+2:]
                        self.bounds[1,2] = float(a)   
                    if "maxk32" in line:
                        ind = line.find('=')                  
                        a = line[ind+2:]
                        self.bounds[1,3] = float(a)   
                    if "maxkr2" in line:
                        ind = line.find('=')                  
                        a = line[ind+2:]
                        self.bounds[1,4] = float(a)   
                    if "maxkr3" in line:
                        ind = line.find('=')                  
                        a = line[ind+2:]
                        self.bounds[1,5] = float(a)   
                    if "maxy" in line:
                        ind = line.find('=')                  
                        a = line[ind+2:]
                        self.bounds[1,6] = float(a)           
                    if "maxb" in line:
                        ind = line.find('=')                  
                        a = line[ind+2:]
                        self.bounds[1,7] = float(a)   
                    if "helpon" in line:
                        ind = line.find('=')
                        a = line[ind+2:]
                        self.helptips.set(int(a))
        except:
            tkMessageBox.showerror(title = 'Error',message='There is no config file present, please save a config file when you want to save your inputs. Reverting to default inputs.')
                            
                
        self.poll()



    def modelpictureframes(self,ModelPictureFrame):
            #set up the model display frame
    
            self.modelcanvas = Canvas(self.ModelPictureFrame,width=404,height=222)
            self.modelcanvas.pack(anchor=NW,padx=2,pady=2)
            self.modelgif = self.modelimages[0]       
            self.modelcanvas.image =self.modelgif
            #define the self call for the image currently on the canvas (to update via swapimage())
            self.image_on_canvas = self.modelcanvas.create_image(0, 0,anchor=NW, image=self.modelgif)
    
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

 

      
def main():
    
    root = Tk() #create GUI app
   
    root.geometry("416x783+10+10") #define GUI shape
    root.iconbitmap('qbio_icon_9nD_icon.ico') #set icon
    #root.resizable(width=False, height=False) #no resizing
    
    app = GUI(root)#def app
    
    root.mainloop()  #do the GUI stuff
    try:    
        root.destroy() #destroy the GUI if the user hits file exit
    except:
        pass # dont if they hit the x button because its already gone


if __name__ == '__main__': #do the main
    main()  