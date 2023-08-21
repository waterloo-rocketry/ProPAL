# NOSThermo : Class to handle nitrous thermophysical properties from the NIST paper (Helmholtz EoS).
# Author : Cristian Bicheru
# NOTE: Cleanup is required. This is a quick proto for Artem to mess around with.
import math
from scipy.optimize import fsolve

Ttp = 182.33
Tbp = 186.68
M = 44.0128
omega = 0.1613

Tc = 309.52
Pc = 7.245
RhocM = 10.27 # mol / L
Rhoc = RhocM * M # kg/m^3​

R = 8.314472

n = [
    0.88045,
    -2.4235,
    0.38237,
    0.068917,
    0.00020367,
    0.13122,
    0.46032,
    -0.0036985,
    -0.23263,
    -0.00042859,
    -0.042810,
    -0.023038
]

def diff(f, x, *K):
    return (f(x+1e-5, *K)-f(x-1e-5, *K))/2e-5

class NOSThermo:
    solverTol = 1e-6

    def __init__(self):
        pass

    def evalHelmPolar(this, delta, tau):
        return n[0] * delta * tau**0.25 + \
            n[1] * delta * tau**1.25 + \
            n[2] * delta * tau**1.5 + \
            n[3] * delta**3 * tau**0.25 + \
            n[4] * delta**7 * tau**0.875 + \
            n[5] * delta * tau**2.375 * math.exp(-delta) + \
            n[6] * delta**2 * tau**2 * math.exp(-delta) + \
            n[7] * delta**5 * tau**2.125 * math.exp(-delta) + \
            n[8] * delta * tau**3.5 * math.exp(-delta**2) + \
            n[9] * delta * tau**6.5 * math.exp(-delta**2) + \
            n[10] * delta**4 * tau**4.75 * math.exp(-delta**2) + \
            n[11] * delta**2 * tau**12.5 * math.exp(-delta**3)
    
    def delHelmDelDelta(this, delta, tau):
        return n[0] * tau**0.25 + \
            n[1] * tau**1.25 + \
            n[2] * tau**1.5 + \
            n[3] * (3 * delta**2) * tau**0.25 + \
            n[4] * (7 * delta**6) * tau**0.875 + \
            n[5] * tau**2.375 * math.exp(-delta) - n[5] * delta * tau**2.375 * math.exp(-delta) + \
            n[6] * (2 * delta) * tau**2 * math.exp(-delta) - n[6] * delta**2 * tau**2 * math.exp(-delta) + \
            n[7] * (5 * delta**4) * tau**2.125 * math.exp(-delta) - n[7] * delta**5 * tau**2.125 * math.exp(-delta) + \
            n[8] * tau**3.5 * math.exp(-delta**2) - n[8] * (2 * delta**2) * tau**3.5 * math.exp(-delta**2) + \
            n[9] * tau**6.5 * math.exp(-delta**2) - n[9] * (2 * delta**2) * tau**6.5 * math.exp(-delta**2) + \
            n[10] * (4 * delta**3) * tau**4.75 * math.exp(-delta**2) - n[10] * (2 * delta**5) * tau**4.75 * math.exp(-delta**2) + \
            n[11] * (2 * delta) * tau**12.5 * math.exp(-delta**3) - n[11] * (3 * delta**4) * tau**12.5 * math.exp(-delta**3)

    def p(this, rho, T):
        delta = rho/Rhoc
        tau = Tc/T

        delHelmdelDelta = this.delHelmDelDelta(delta, tau)

        return rho * R * T * (1 + delta*delHelmdelDelta)/M * 1000


    # fsolve is slow and doesnt work well here
    def rho(this, p, T):
        return fsolve(lambda x: this.p(x, T) - p, 1000)[0]


    def Z(this, rho, T):
        delta = rho/Rhoc
        tau = Tc/T

        return 1 + delta*this.delHelmDelDelta(delta, tau)