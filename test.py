# Source document
# http://www.aspirespace.org.uk/downloads/Modelling%20the%20nitrous%20run%20tank%20emptying.pdf

from contextlib import suppress
from math import sqrt
from math import exp, pi
import matplotlib.pyplot as plt
from numpy import interp

# This doesn't quite work yet, will figure out why later
from NOSThermo import NOSThermo

# This also errors out when I try it
from tank_model import z_factor


# Project imports



DEBUG_VERBOSE = False

tCrit = 309.57
rhoCrit = 452.0
pCrit = 72.51
ZCrit = 0.28
# NOS ratio os specific heats
ratio_of_specific_heats_gamma = 1.3



class NOS:

    def sus_Z(self, P):
        """
        The incredibly suspicious way that Rick Newlands chooses to calculate his Z, but it is super simple and
        computationally forgiving, is being used until we optimize stuff. Also nothing else works super well yet
        """
        return interp(P, [0.0, pCrit], [1.0, ZCrit])

    # Nitrous oxide vapour pressure, Bar
    def NOS_vapor_pressure(self, T_Kelvin):
        """
        Shamelessly stolen from the Rick Newlands paper (converted from cpp to python 
        with some var names cleaned up)
        """

        p = [1.0, 1.5, 2.5, 5.0]
        b = [-6.71893, 1.35966, -1.3779, -4.051]

        Tr = T_Kelvin / tCrit;
        Tr_minus_one = 1.0 - Tr
        acc = 0.0

        for idx in range(4):
            acc += b[idx] * pow(Tr_minus_one, p[idx])

        return pCrit * exp((acc / Tr))

    # Nitrous oxide saturated liquid density, kg/m3 
    def NOS_liquid_density_saturated(self, T_Kelvin):
        """
        Shamelessly stolen from the Rick Newlands paper (converted from cpp to python 
        with some var names cleaned up)
        """
    
        b = [1.72328, -0.8395, 0.5106, -0.10412]

        Tr = T_Kelvin / tCrit
        Tr_minus_one = 1.0 - Tr
        acc = 0.0
        
        for idx in range(4):
            acc += b[idx] * pow(Tr_minus_one,((idx+1) / 3.0))

        return rhoCrit * exp(acc)

    def NOS_vapor_density_saturated(self, T_Kelvin):        
        """
        Shamelessly stolen from the Rick Newlands paper (converted from cpp to python 
        with some var names cleaned up)
        """

        b = [-1.009, -6.28792, 7.50332, -7.90463, 0.629427]

        Tr = T_Kelvin / tCrit
        Tr_adjusted = (1.0 / Tr) - 1.0
        acc = 0.0
        for idx in range(5):
            acc += b[idx] * pow(Tr_adjusted,((idx+1) / 3.0))
        
        return rhoCrit * exp(acc)
        
    def NOS_enthlapy_of_vaporization(self, T_Kelvin):
        """
        Shamelessly stolen from the Rick Newlands paper (converted from cpp to python 
        with some var names cleaned up)
        """

        bL = [-200.0, 116.043, -917.225, 794.779, -589.587]
        bV = [-200.0, 440.055, -459.701, 434.081, -485.338]

        Tr = T_Kelvin / tCrit
        rab = 1.0 - Tr
        acc_L = bL[0]
        acc_V = bV[0]

        #for (dd = 1; dd < 5; dd++)


        for idx in range(1,5):
            acc_L += bL[idx] * pow(rab,(idx / 3.0)); # saturated liquid enthalpy 
            acc_V += bV[idx] * pow(rab,(idx / 3.0)); # saturated vapour enthalpy 
        

        bob = (acc_V - acc_L) * 1000.0; # net during change from liquid to vapour 
        return(bob)
    
    def NOS_isobaric_heat_capacity(self, T_Kelvin):
        """
        Shamelessly stolen from the Rick Newlands paper (converted from cpp to python 
        with some var names cleaned up)
        """
        b = [2.49973, 0.023454, -3.80136, 13.0945, -14.518]
        Tr = T_Kelvin / tCrit
        Tr_minus_one = 1.0 - Tr
        acc = 1.0 + b[1] / Tr_minus_one

        #for (dd = 1; dd < 4; dd++)
        for idx in range(1,4):
            acc += b[(idx+1)] * pow(Tr_minus_one,idx)

        result = b[0] * acc * 1000.0 # convert from KJ to J */
        return result

    

    # def basic_compressibility(T, P):
    #     return 1

    def calculate_CC_pressure(self, tank_pressure):
        """
        Temporary ballpark estimate until we get an actual Combustion Chamber model running
        """
        return tank_pressure*0.5
    
    def calculate_loss_factor(self, loss_coeff):
        """
        Calculate the flow losses 
        """

        NUM_ORIFICES = 36
        ORIFICE_DIAM = 0.003 # 3mm

        # not currently used, we started picking values to match massflow, see comment above func. call
        MAGIC_ASPIRESPACE_LOSS_COEFF = 2*30 

        # Diameter is not squared, this may be wrong, needs further investigation

        # ^ That was actually corrected, but all of this is still incredibly sketch, needs more
        # investigation

        loss_factor = loss_coeff/((0.25*pi*(ORIFICE_DIAM**(2))*NUM_ORIFICES)**2)
        return loss_factor

    def calc_densities(self):
        self.density_liquid = self.NOS_liquid_density_saturated(self.temperature)
        self.density_vapor = self.NOS_vapor_density_saturated(self.temperature)

    def calc_masses(self):
        vapor_volume = self.ullage * self.volume
        liquid_volume = self.volume  - vapor_volume

        self.mass_liquid = liquid_volume * self.density_liquid
        self.mass_vapor = vapor_volume * self.density_vapor

    def calculate_injector_massflow(self, P_tank, P_comb_chamber, target_density, loss_coeff):
        pressure_drop = P_tank - P_comb_chamber

        mass_flowrate = sqrt(2 * target_density * pressure_drop*100000 / self.calculate_loss_factor(loss_coeff))
        return mass_flowrate
    

    def execute_repressurization_cooling(self, mass_vaporized):
        """
        Account for the energy required to vaporize the gas as the tank empties in the saturated phase
        """
        curr_enthlapy = self.NOS_enthlapy_of_vaporization(self.temperature)
        curr_C_liquid = self.NOS_isobaric_heat_capacity(self.temperature)
        
        # print('Mass Vaporized: ' + str(mass_vaporized))

        curr_delta_Q = mass_vaporized*curr_enthlapy 
        curr_delta_T =  - (curr_delta_Q) / (self.mass_liquid * curr_C_liquid)



        self.temperature += curr_delta_T
        self.pressure = self.NOS_vapor_pressure(self.temperature)
        self.calc_densities()

    
    def calc_temperature_during_vap_only(self, T_prev, m_prev, Z_prev, Z, m):
        """
        Use the isentropic gas depressurization formulas to calculate temperature
        """

        eqn_exponent = (ratio_of_specific_heats_gamma - 1)
        return T_prev*(((Z*m)/(Z_prev*m_prev))**(eqn_exponent))
    
    def calc_pressure_during_vap_only(self, T_prev, P_prev, T):
        """
        Use the isentropic gas depressurization formulas to calculate pressure
        """

        eqn_exponent = (ratio_of_specific_heats_gamma - 1)/ratio_of_specific_heats_gamma
        return ((T)/(T_prev))**(1/(eqn_exponent))*P_prev

    def calc_density_during_vap_only(self, T_prev, rho_prev, T):
        """
        Use the isentropic gas depressurization formulas to calculate vapor density
        """

        eqn_exponent = 1.0/(ratio_of_specific_heats_gamma - 1.0)
        return rho_prev * ((T/T_prev)**(eqn_exponent))


    def execute_vapor_phase_state(self, suppress_prints=True):
        """
        Handles the ideal isentropic venting state, which requires an iterative process for finding the right Z.

        This process is detailed throughly by Rick Newlands, but it may be worth improving in the future.
        """


        previous_Z = self.sus_Z(self.pressure)
        
        #thermo = NOSThermo()
        #previous_Z = thermo.Z(self.density_vapor, self.temperature)

        previous_vapor_mass = self.mass_vapor + self.interval_delta_mass
        guess_Z = previous_Z
        exit_flag = False

        z_iter_count = 0
        conv_step = 10/9
        ITER_LIMIT = 100000

        aim = 0

        while not exit_flag:
            iter_T = self.calc_temperature_during_vap_only(self.temperature, previous_vapor_mass,
                                                    previous_Z, guess_Z,
                                                    self.mass_vapor) 

            iter_P = self.calc_pressure_during_vap_only(self.temperature, self.pressure, iter_T)
            
            # iter_Z = z_factor(iter_T, iter_P, perform_plotting=False)
            iter_Z = self.sus_Z(iter_P)
            # iter_Z = thermo.Z(self.density_vapor, iter_T)

            if (z_iter_count  % 5 == 0 and not suppress_prints):
                print('Iteration T:' + str(iter_T) + ' Iteration P:' +\
                        str(iter_P) + ' delta_Z: ' + str(guess_Z - iter_Z) + ' Iter Count: ' + str(z_iter_count))

            # Convergence achieved
            if abs(iter_Z - guess_Z) < 0.000001:
                # print('Z Guess: ' + str(iter_Z))
                exit_flag = True
                self.pressure = iter_P
                self.previous_temperature = self.temperature
                self.temperature = iter_T
                # print("Number of Z iterations: " + str(z_iter_count) + " Calculated Z value: " + str(iter_Z))
                return True

            # Adjust guess based on relative size
            oldAim = aim
            if iter_Z > guess_Z:
                guess_Z = guess_Z*conv_step
                aim = -1
            elif iter_Z < guess_Z:
                guess_Z = guess_Z/conv_step 
                aim = 1
            
            if aim == -oldAim:
                conv_step = sqrt(conv_step)


            # Check for convergence failure
            z_iter_count += 1
            if z_iter_count > ITER_LIMIT:
                print('Failed to converge while iterating Z')
                exit_flag = True
                return False




    def execute_vapack(self, time_step, suppress_prints=False):
        """
        Execute a single cycle of tank blowdown under self-pressurization. 

        This is referred by Rick Newlands as 'Vapak'
        """

        curr_cc_pressure = self.calculate_CC_pressure(self.pressure)

        # These numbers were adjusted ad-hoc while we were in the bar with help from Cristian B. and Aaron L. to roughly
        # match the massflow of our engine, this will need to be looked at later in more detail
        if self.mass_liquid > 0:
            interval_massflow = self.calculate_injector_massflow(self.pressure, curr_cc_pressure, self.density_liquid, 80)
        else:
            interval_massflow = self.calculate_injector_massflow(self.pressure, curr_cc_pressure, self.density_vapor, 10)

        # As per the paper, using Adams-Bashforth 2nd Order Equation
        # Source: https://en.wikipedia.org/wiki/Linear_multistep_method

        previous_massflow = self.massflow
        previous_delta_mass = self.interval_delta_mass

        # Aspirespace paper doesn't use the first term of the Adams-Bashford, and I can't get it to work like that anyway
        # This is now essentially copypasta from the example code (link restated for convenience) 
        # http://www.aspirespace.org.uk/downloads/Modelling%20the%20nitrous%20run%20tank%20emptying.pdf
        # (Page 14)

        delta_system_mass = time_step * ((3/2)*interval_massflow - (1/2)*previous_massflow)

        # Alternative (non-integrated) form also seems like it might work, leaving it here:
        #delta_system_mass = time_step * interval_massflow

        self.massflow = interval_massflow
        self.interval_delta_mass = delta_system_mass

        if self.mass_liquid <= 0: # Final stage of burn; vapor-only expansion
            
            
            self.mass_liquid = 0
            self.mass_vapor -= self.interval_delta_mass
            no_error = self.execute_vapor_phase_state(suppress_prints=True)
            self.density_vapor = self.calc_density_during_vap_only(self.previous_temperature, self.density_vapor, self.temperature)

            if self.pressure <= 2.5: # Random bullshit number pulled out of thin air by Rick Newlands to stop the burn at 
                return 2

            if not no_error: # There be errors
                return -3


        else: # Standard first burn stage, gas-liquid equilibrium in tank

            # Prior to accounting for nitrous vaporization to account 
            self.prev_mass_liquid = self.mass_liquid
            self.mass_liquid -= delta_system_mass 
            curr_system_mass = self.mass_liquid + self.mass_vapor

            if not suppress_prints:
                print('Current system mass: ' + str(curr_system_mass), end = ' ')
                print('Calculated Vapor Density: ' + str(self.density_vapor), end = ' ')
                print('Calculated Liquid Density ' + str(self.density_liquid))

            # Actual mass of the liquid in the tank determined based on liquid-gas equilibrium in the tank
            # given reduced overall system mass    
            true_liquid_mass = (self.volume  - (curr_system_mass / self.density_vapor))/ \
                    (self.density_liquid**(-1) - self.density_vapor**(-1))

            if not suppress_prints:
                print('Prevaporized liquid mass: ' + str(self.mass_liquid), end = ' ')
                print('True liquid mass: '+ str(true_liquid_mass))

            # The difference in masses is the liquid that would be vporized to maintain vapor pressure
            mass_vaporized = self.mass_liquid - true_liquid_mass


            self.mass_vapor +=  mass_vaporized
            self.mass_liquid = true_liquid_mass

            # Record first-order lagging vapor mass value to account for vaporization time of the nitrous
            # In the same way as is done in the aspirespace paper
            # Also without this the numerical model eats shit and begins to oscillate and die

            self.lagging_vaporized_mass += (time_step/0.15) * (mass_vaporized - self.lagging_vaporized_mass)
            self.execute_repressurization_cooling(self.lagging_vaporized_mass)


            if self.mass_liquid <= 0:
                print("Liquid Exhausted")
                self.mass_liquid = 0

        if not suppress_prints:
            print("Iter CC Pressure: " + str(curr_cc_pressure), end = ' ')
            print('Massflow: ' + str(self.massflow))
            print("Iter Fluid Temperature: " + str(self.temperature))
            
        
        if self.mass_vapor <= 0:
            return -1
        elif curr_cc_pressure < 1: 
            return 0 # CC pressure below atmospheric, gas expansion ceases

        return 1 # Ready to continue iterating


    def __init__(self, V_0, P_0, T_0, m_0, basic_ullage) -> None:
        """
        Sets the intial conditions. 

        NOTE/TODO: the pressure and mass are a function of the other values, and they are not actually necessary.
        However, they might still have utility for providing values/performing checks. Decide this later. 
        """

        self.ullage = basic_ullage
        self.temperature = T_0
        self.volume = V_0
        self.pressure = self.NOS_vapor_pressure(T_0)
        
        self.density_liquid = -1
        self.density_vapor = -1

        self.calc_densities()

        self.mass_liquid = -1
        self.mass_vapor = -1

        self.calc_masses()

        self.interval_delta_mass = 0
        self.massflow = 0
        self.lagging_vaporized_mass = 0






if __name__ == "__main__":
    exit_flag = False
    iter_count = 0

    # Initialize NOS model
    NOS_model = NOS(0.04, 55, 288, 40, 0.15)

    sim_time = 0
    timestamps = []
    pressure_probe = []
    temp_probe = []
    mass_flow_probe = []
    mass_probe = []
    liquid_mass_probe = []
    vapor_mass_probe = []
    lagging_vaporized_probe = []

    TIME_STEP = 0.0001
    ITER_LIMIT = 1000000

    while (iter_count < ITER_LIMIT and not exit_flag):
        suppress_prints = True
        if iter_count % 1000 == 0:
            suppress_prints = False
            print('Simulation time: ' + str(sim_time), end=' | ')

        status = NOS_model.execute_vapack(TIME_STEP,suppress_prints=suppress_prints)

        if status <= 0 or status == 2:
            exit_flag = True

        sim_time += TIME_STEP
        sim_time = round(sim_time, 7)

        timestamps.append(sim_time)
        pressure_probe.append(NOS_model.pressure)
        temp_probe.append(NOS_model.temperature)
        mass_flow_probe.append(NOS_model.interval_delta_mass*100) # To scale simulation better
        mass_probe.append(NOS_model.mass_liquid + NOS_model.mass_vapor)
        liquid_mass_probe.append(NOS_model.mass_liquid)
        vapor_mass_probe.append(NOS_model.mass_vapor)
        lagging_vaporized_probe.append(NOS_model.lagging_vaporized_mass)


        iter_count += 1

        if exit_flag:
            print("Terminated ", end = '')
            if status == -3:
                print('due to convergence failure while calculating Z')
            elif status == 0:    
                print("due to standard exit condition (CC pressure below atmospheric)") 
            elif status == 2:
                print("due to Rick Newland's Magic 2.5 Bar vapor pressure criteria")   
            else:
                print('due to nonstandard exit condition')

            print("Simulation time: " + str(sim_time))

            if (DEBUG_VERBOSE):
                print(pressure_probe)
                print(temp_probe)
                print(mass_flow_probe) 

            plt.plot(timestamps, pressure_probe, label = 'Pressure')
            plt.plot(timestamps, temp_probe, label = 'Temperature')
            plt.plot(timestamps, mass_flow_probe, label = 'Mass Flow')
            plt.plot(timestamps, mass_probe, label = 'Mass')
            plt.plot(timestamps, liquid_mass_probe, label = 'Liquid Mass')
            plt.plot(timestamps, vapor_mass_probe, label = 'Vapor Mass')
            plt.plot(timestamps, lagging_vaporized_probe, label = 'Lagging Vapor Mass')
            

            plt.legend()

            plt.show()

            
    if iter_count >= (ITER_LIMIT - 1):
        print("Main loop iteration count exceeded")





        

