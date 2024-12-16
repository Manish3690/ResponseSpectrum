import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from nonlinearanalysis import nonLinearAnalysis
import math





N=1
KN=1e3*N
m=1
cm=1.e-2*m
mm=1e-3*m
Pa=N/m**2
MPa=1e6*Pa



E= 2.e5*MPa
I= 301*cm**4
di=102*mm
do=115*mm
L=0.1*m
M= 1650.0 #kg  not used the mass of the column itself as mass is not much while compared to lumped mass
Tn=0
Time_array=np.zeros(600)
u_array=np.zeros(600)
a_array=np.zeros(600)

damp_ratio= [0.05, 0.1, 0.15,0.2]



earthquake = pd.read_csv("C:/Users/acer/Desktop/ResponseSpectra/gorkha_dt 0.005.txt")
p_df= earthquake*M*9.81
p=list(p_df.values.flatten())
p=np.array(p)
data_spectrum= pd.DataFrame()

for j in range(len(damp_ratio)):
    L=0.1*m
    for i in range(600):
        non= nonLinearAnalysis(p, M, damp_ratio[j], E, I, L, 0.005)
        u_array[i],Time_array[i] = non.Newmark_normal()
        a_array[i]= ((2*math.pi/(Time_array[i]))**2)*u_array[i]
        L+=0.05
        print(Time_array[i])
        
        
    data_spectrum["Time"]=Time_array
    data_spectrum[f"{damp_ratio[j]}"]= a_array
       

#%%

data_spectrum.to_csv("data_spectrum.csv")
    
plt.plot(Time_array, data_spectrum.iloc[:,1], color='red')
plt.plot(Time_array, data_spectrum.iloc[:,2], color='green')
plt.plot(Time_array, data_spectrum.iloc[:,3], color='blue')
plt.plot(Time_array, data_spectrum.iloc[:,4], color='yellow')
# plt.legend()
# plt.axis([0.,1.,0.,6.])
plt.show()
plt.savefig("SampleGorkhaPseudoAccSpectrum.jpg")


print(max(Time_array))