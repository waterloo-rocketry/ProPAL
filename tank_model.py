from scipy.optimize import fsolve
import numpy as np
from scipy.interpolate import make_interp_spline
import matplotlib.pyplot as plt

R = 0.08206 #universal constant 
Pc = 72.420429 #unit in atm
Tc = 309.56 #unit in kelvin 

a = 27 * R**2 * Tc**2 / (Pc * 64)
b = R * Tc / (8 * Pc)



def z_factor(my_temp,my_pressure, perform_plotting=True):

    Pr = my_pressure/Pc #reduced pressure 
    Tr= my_temp/Tc #reduced temp
    
    def function_gas(v,Pr):
        p = Pr * Pc
       
        return p * (v - b) - (R * my_temp)  +  a / v**2 * (v - b) #non idead gas equation 

    Pr_range = np.linspace(0.05, 7) #defining the range of the solution 
    v = [fsolve(function_gas,1, args=(Pr,))[0] for Pr in Pr_range] # solve for the Volume in non-ideal gas formula for all the pressure range 
    P_range = Pr_range * Pc #pressure range to find all the z factors with specifed range
    Z = P_range * v/ (R * my_temp) #calculate Z factor

    if not perform_plotting:
        return Pr_range[0]*Pc*v[0]/(R * my_temp)

    #print(Pr_range)
    ax = plt.subplot(1, 1, 1)

# Major ticks every 1 unit, minor ticks every 0.5 for x axis 
    major_ticks = np.arange(0, 8, 1)
    minor_ticks = np.arange(0, 8, 0.5)
    major_ticks2 = np.arange(0, 2, 0.2)
    minor_ticks2 = np.arange(0, 2, 0.125)
    ax.set_xticks(major_ticks)
    ax.set_xticks(minor_ticks, minor=True)
    ax.set_yticks(major_ticks2)
    ax.set_yticks(minor_ticks2, minor=True)

    # And a corresponding grid
    ax.grid(which='both')

    # Or if you want different settings for the grids:
    ax.grid(which='minor', alpha=0.2)
    ax.grid(which='major', alpha=0.5)
    #curve smoother 
    B_spline_coeff = make_interp_spline(Pr_range, Z)
    X_Final = np.linspace(Pr_range.min(), Pr_range.max(), 250)
    Y_Final = B_spline_coeff(X_Final)
    plt.plot(X_Final,Y_Final)
    plt.legend()
    plt.xlabel('$P_r$ atm')
    plt.ylabel('Z')
    plt.xlim([0, 7])
    plt.ylim([0, 1.2])

    #print(P_range[0]*v[0]/(R * my_temp))
    return Pr_range[0]*Pc*v[0]/(R * my_temp)
    

if __name__ == '__main__':
    for t in range(330,650,20):
        z_factor(t,1)
    
    plt.show()




