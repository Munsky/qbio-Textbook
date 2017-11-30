
import numpy as np
import scipy as sci
from scipy.integrate import odeint 
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