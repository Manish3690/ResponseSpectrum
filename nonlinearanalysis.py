import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os




class nonLinearAnalysis:
    def __init__(self,p,mass,damp_ratio,E,I,L,del_t):
        self.p = p
        self.mass = mass
        self.damp_ratio = damp_ratio
        self.del_t = del_t
        
        self.E=E
        self.I = I
        self.L =L
        
        self.k=3*self.E*self.I/(self.L**3)
        self.nat_frrq=np.sqrt(self.k/self.mass)
        self.Tn=2*np.pi/self.nat_frrq
        self.c = 2*self.mass*self.nat_frrq* self.damp_ratio
        pass
    

    

    
    def Newmark_normal(self):
        """
        Perform dynamic analysis using the Newmark-beta method.
    
        Parameters:
        earthquake (list or np.array): Ground acceleration time history.
        mass (float): Mass of the system.
        damping_ratio (float): Damping ratio.
        nat_frq (float): Natural frequency of the system.
        k (float): Stiffness of the system.
        del_t (float): Time step.
    
        Returns:
        u_max (float): maximum displacement.
        """
        # Newmark parameters for average acceleration
        alpha = 0.5
        beta = 0.25

    
    
        #a1,a2,a3
    
        
        a1 = self.mass * 4 / (self.del_t)**2 + self.c * 2 / self.del_t
        a2 = self.c + 4 * self.mass / self.del_t
        a3 = self.mass

        k_app = self.k + a1

        n_steps = len(self.p)
        u = np.zeros(n_steps + 1)  # Displacement
        vel = np.zeros(n_steps + 1)  # Velocity
        acc = np.zeros(n_steps + 1)  # Acceleration
        p1 = np.zeros(n_steps + 1)

        for i in range(n_steps):
            acc[i + 1] = (self.p[i] - self.c * vel[i] - self.k * u[i]) / self.mass
            p1[i + 1] = self.p[i] + a1 * u[i] + a2 * vel[i] + a3 * acc[i]
            u[i + 1] = p1[i + 1] / k_app
            vel[i + 1] = 2 * (u[i + 1] - u[i]) / self.del_t - vel[i]
            acc[i + 1] = 4 * (u[i + 1] - u[i]) / self.del_t**2 - 4 * vel[i] / self.del_t - acc[i]

        # return max(abs(u)), self.Tn
        return u, vel, acc

    
    def nonLinear(self,fy):
        
        u=[0]*len(self.p)
        vel=[0]*len(self.p)
        acc=[0]*len(self.p)
        fs=[0]*len(self.p)
        k_T=[0]*len(self.p)
        k_T[0]=self.k
        p_1=[0]*len(self.p)
        k_T1=np.zeros(len(self.p))
        
        u_linear, vel, acc = self.Newmark_normal()
        
        a1= 4*self.mass/(self.del_t)**2 + 2*self.c/self.del_t
        a2= 4*self.mass/self.del_t + self.c
        a3= self.mass
        
        
        for i in range(len(self.p)-1):
          u[i+1]=u_linear[i]
          fs[i+1]=fs[i]
          k_T[i+1]= k_T[i]
        
        
          p_1[i+1]= self.p[i+1] + a1* u[i] + a2 * vel[i] + a3 * acc[i]
          Residue = p_1[i+1] - fs[i+1] - a1 * u[i+1]
          j=0
          # print(p_1[i+1], Residue)
          while abs(Residue) >= 1.e-5:
            j+=1
            # print(j)
            k_T1[i+1] = k_T[i+1] + a1
            del_u = Residue/ k_T1[i+1]
            # print(del_u)
            u[i+1] = u[i+1]+del_u
            fs[i+1] = fs[i] + (u[i+1] - u[i])* self.k
            if abs(fs[i+1])  >= fy:
              fs[i+1] = 30*np.sign(fs[i+1])
              k_T[i+1]=0
            Residue= p_1[i+1] - fs[i+1] - a1 * u[i+1]
            print(j, Residue, u[i+1], fs[i+1])
            
        return u,fs
            
            
            
    
    
    
    
    
    
    
    
