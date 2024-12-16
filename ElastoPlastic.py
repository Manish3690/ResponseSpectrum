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
L=3*m
M= 1650.0



df= pd.read_csv('gorkha_dt 0.005.txt')
p_df=df*m*981
p=list(p_df.values.flatten())
p=np.array(p)
print(p)



non= nonLinearAnalysis(p, M, 0.05, E, I, L, 0.005)
u,fs= non.nonLinear(30)
plt.plot(u,fs)
plt.show()



