from scipy.optimize import fsolve
import numpy as np
import matplotlib.pyplot as plt

R = 0.08206 #universal constant 
Pc = 72.420429 #unit in atm
Tc = 309.56 #unit in kelvin 

a = 27 * R**2 * Tc**2 / (Pc * 64)
b = R * Tc / (8 * Pc)



def z_factor(my_temp,my_pressure):

    Pr = my_pressure/Pc #reduced pressure 
    print(Pr)
    Tr= my_temp/Tc #reduced temp
    
    def function_gas(V,Pr):
        P = Pr * Pc
       
        return P * (V - b) - (R * my_temp)  +  a / V**2 * (V - b) #non idead gas equation 

    Pr_range = np.linspace(0.1, 10) #defining the range of the solution 
    V = [fsolve(function_gas,4, args=(Pr,))[0] for Pr in Pr_range] # solve for the Volume in non-ideal gas formula for all the pressure range 

    
    P_range = Pr_range * Pc #pressure range to find all the z factors with specifed range
    Z = P_range * V / (R * my_temp) #calculate Z factor 

    plt.plot(Pr_range, Z)
    plt.xlabel('$P_r$ atm')
    plt.ylabel('Z')
    plt.xlim([0, 10])
    plt.ylim([0, 2])
