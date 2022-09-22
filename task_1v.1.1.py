
#import required library 
from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt
import pint 


#defining the constants 

#molar mass
M = 44.028 #unit in mole 
#critical temp 
Tc = 309.56 #unit in kelvin
#critical pressure 
Pc = 73.38 # unit in bar
#universal gas constatn 
Ru= 8.3144598 #unit in J.mol.k
#specific gas constant 
R = Ru/M

Temp = 90 # what is going to be passed on as temp 
Pressure = 20 #pressure is going to be passed as input 
#find reduced temp and pressure 
Tr = Temp/Tc
Pr = Pressure/Pc




a = 27 * R**2 * Tc**2 / (Pc * 64)
b = R * Tc / (8 * Pc)


def objective(V, Pr):
    P = Pr * Pc
    T = Tr * Tc
    return P * (V - b) - (R * T)  +  a / V**2 * (V - b)


Pr_range = np.linspace(0.1, 10)
V = [fsolve(objective, 3, args=(Pr,))[0] for Pr in Pr_range]

T = Tr * Tc
P_range = Pr_range * Pc
Z = P_range * V / (R * T)

plt.plot(Pr_range, Z)
plt.xlabel('$P_r$')
plt.ylabel('Z')
plt.xlim([0, 10])
plt.ylim([0, 2])
