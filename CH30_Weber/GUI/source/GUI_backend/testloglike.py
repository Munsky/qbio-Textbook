# -*- coding: utf-8 -*-
"""
Created on Mon Jan 30 13:31:57 2017

@author: William
"""
try:
    import GUI_backend.Gene_Model_Class as Gene_Model_Class #backend class 
except:
    import Gene_Model_Class
import numpy as np

GeneModel = Gene_Model_Class.GeneModel()

N = 20
states = np.linspace(0,N,N)
x0 = np.zeros((1,3*N)).T
x0[0] =1
x0[1] = 0
x0[2] =0
x0[3] = 0
t = np.linspace(0,100,100)
Parameters = np.array([1,  1,   .5,   .5, 2, .5,  .2,1])
mRNA_Hist = GeneModel.Import_mRNA("C:/Users/William/Desktop/mRNA_Hist.csv")
soln,N = GeneModel.FSPCalc_ODE_prealloc(Parameters,'',0,N,x0,t,.1)


    #reshape our for easier computing              
logsum = 0

for i in range(0,len(mRNA_Hist[0])):
    
    nhist = np.histogram(mRNA_Hist.T[i],bins = np.linspace(min(mRNA_Hist.T[i]),max(mRNA_Hist.T[i]),(max(mRNA_Hist.T[i])-min(mRNA_Hist.T[i])+1)))
    x,RNAprob = GeneModel.FSP_Return_RNA_Prob(soln[i],N)
   
    for j in range(0,len(nhist[0])):

        if j > len(RNAprob)-1:
            RNAprob = np.append(RNAprob,1e-10)
            

        if RNAprob[j] <=0:
            RNAprob[j] = 1e-10            
    
        logsum = logsum + nhist[0][j]*np.log(RNAprob[j])
    
        
    
    
print -logsum



Bounds = GeneModel.Get_Par_Bounds(2)
parameters = np.array([1,  1,   .5,   .5, 2, .5,  .2,1])
NBurn = 10  
NChain = 100
NThin = 10
MutRate = .2

x_in = [1,0,0,0]
TimeVar = GeneModel.TimeVar_Create('(1-cos(2*3.14/30*t))*(t>5)*(t<70)')
Data_Set = GeneModel.Import_mRNA('C:/Users/William/Desktop/mRNA_Hist.csv')

tol = .1
Function = lambda x,N: FSP_fit_fun(x,2,t,x_in,Data_Set,TimeVar,tol,N)

'''

if f(x+i) < f(xi)
Ncurv = Nexpand
maccept = maccept + 1

if maccept = 10
 Ncur = Ninit
'''


def FSP_fit_fun(Parameters,case,t,x_in,mRNA_Hist,TimeVar,tol,N):
    x0 = np.zeros((1,3*N)).T
    x0[0] = x_in[0]
    x0[1] = x_in[1]
    x0[2] = x_in[2]
    x0[3] = x_in[3]
    
    
    soln,N = GeneModel.FSPCalc_ODE_prealloc(Parameters,TimeVar,case,N,x0,t,tol)
    
  
    logsum = 0

    for i in range(0,len(mRNA_Hist[0])):
        nhist = np.histogram(mRNA_Hist.T[i],bins = np.linspace(min(mRNA_Hist.T[i]),max(mRNA_Hist.T[i]),(max(mRNA_Hist.T[i])-min(mRNA_Hist.T[i])+1)))
      
        x,RNAprob = GeneModel.FSP_Return_RNA_Prob(soln[i],N)
        for j in range(0,len(nhist[0])):

            if j > len(RNAprob)-1:
                RNAprob = np.append(RNAprob,1e-10)
                
    
            if RNAprob[j] <=0:
                RNAprob[j] = 1e-10       
                        
             
                
            logsum = logsum + nhist[0][j]*np.log(RNAprob[j])
                
    return -logsum,N


def MetroplisHastings_FSP(Function,Init_Parameters,Bounds,N_Burn,N_Chain,N_Thin,Mut_rate,N):
    x0 = np.copy(Init_Parameters)
    x = np.copy(x0)#self.make_binaryvec(x0,Bounds[0,:],Bounds[1,:],N,M) #convert the intial parameter guess to binary
    xbest = np.copy(x) #intialize the xbest storage variable
    f = Function #define the function
    #xfloat = self.make_floatvec(x,Bounds[0,:],Bounds[1,:],N,M) #intialize xfloat
    fx,N = f(x,N) #intialize fx
    fbest = np.copy(fx) #initalize fbest
    funchain = np.array([]) #intialize output variables
    parchain = Init_Parameters
    Ncounter = 0
    Nlist = np.array([]) 
    Ncur = N
    for i in range(-N_Burn,N_Chain): #for the length of Nchain, record
        
        for j in range(0,N_Thin):

            xnew = np.copy(x) #copy the binary vector
            for k in range(0,len(xnew)):
                if np.random.random() < Mut_rate:
                    xnew[k] = (Bounds[1,k] - Bounds[0,k])*np.random.random()+Bounds[0,k]

                   
            Nold = N
            print N
            fnew,N = f(xnew,N)
            if N == Nold:
                Ncounter += 1
                if Ncounter > 4:
                    N = int(np.ceil(.7*N))

        
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
                    Nlist = np.append(Nlist,N)
                    
                    
    return funchain, parchain, fbest, xbest, Nlist                      

Best_Par_MetHaste = MetroplisHastings_FSP(Function,parameters,Bounds,NBurn,NChain,NThin,MutRate,20)

