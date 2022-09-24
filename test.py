# Source document
# http://www.aspirespace.org.uk/downloads/Modelling%20the%20nitrous%20run%20tank%20emptying.pdf

from math import sqrt
from math import exp, pi
import matplotlib.pyplot as plt

DEBUG_VERBOSE = False

tCrit = 309.57
rhoCrit = 452.0
pCrit = 72.51

# NOS ratio os specific heats
ratio_of_specific_heats_gamma = 1.3



class NOS:


    # Nitrous oxide vapour pressure, Bar
    def NOS_vapor_pressure(self, T_Kelvin):

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
    
        b = [1.72328, -0.8395, 0.5106, -0.10412]

        Tr = T_Kelvin / tCrit
        Tr_minus_one = 1.0 - Tr
        acc = 0.0
        
        for idx in range(4):
            acc += b[idx] * pow(Tr_minus_one,((idx+1) / 3.0))

        return rhoCrit * exp(acc)

    def NOS_vapor_density_saturated(self, T_Kelvin):        
        b = [-1.009, -6.28792, 7.50332, -7.90463, 0.629427]

        Tr = T_Kelvin / tCrit
        Tr_adjusted = (1.0 / Tr) - 1.0
        acc = 0.0
        for idx in range(5):
            acc += b[idx] * pow(Tr_adjusted,((idx+1) / 3.0))
        
        return rhoCrit * exp(acc)
        
    def NOS_enthlapy_of_vaporization(self, T_Kelvin):

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
        b = [2.49973, 0.023454, -3.80136, 13.0945, -14.518]
        Tr = T_Kelvin / tCrit
        Tr_minus_one = 1.0 - Tr
        acc = 1.0 + b[1] / Tr_minus_one

        #for (dd = 1; dd < 4; dd++)
        for idx in range(1,4):
            acc += b[(idx+1)] * pow(Tr_minus_one,idx)

        result = b[0] * acc * 1000.0 # convert from KJ to J */
        return result

    

    def basic_compressibility(T, P):
        return 1

    def calculate_CC_pressure(self, tank_pressure):
        return tank_pressure*0.5
    
    def calculate_loss_factor(self):
        NUM_ORIFICES = 36
        ORIFICE_DIAM = 0.003 # 3mm
        MAGIC_ASPIRESPACE_LOSS_COEFF = 2

        loss_factor = MAGIC_ASPIRESPACE_LOSS_COEFF/((0.25*pi*ORIFICE_DIAM*NUM_ORIFICES)**2)
        return loss_factor

    def calc_densities(self):
        self.density_liquid = self.NOS_liquid_density_saturated(self.temperature)
        self.density_vapor = self.NOS_vapor_density_saturated(self.temperature)

    def calc_masses(self):
        vapor_volume = self.ullage * self.volume
        liquid_volume = self.volume  - vapor_volume

        self.mass_liquid = liquid_volume * self.density_liquid
        self.mass_vapor = vapor_volume * self.density_vapor

    def calculate_injector_massflow(self, P_tank, P_comb_chamber):
        pressure_drop = P_tank - P_comb_chamber

        mass_flowrate = sqrt(2 * self.density_liquid * pressure_drop / self.calculate_loss_factor())
        return mass_flowrate
    

    def execute_repressurization_cooling(self, mass_vaporized):
        curr_enthlapy = self.NOS_enthlapy_of_vaporization(self.temperature)
        curr_C_liquid = self.NOS_isobaric_heat_capacity(self.temperature)
        
        curr_delta_Q = mass_vaporized*curr_enthlapy 
        curr_delta_T =  - (curr_delta_Q) / (self.mass_liquid * curr_C_liquid)

        self.temperature += curr_delta_T
        self.pressure = self.NOS_vapor_pressure(self.temperature)
        self.calc_densities()

    def execute_vapor_phase_state(self):
        initial_z = self.basic_compressibility(self.temperature, self.pressure)

        return True


    def execute_vapack(self, time_step):

        curr_cc_pressure = self.calculate_CC_pressure(self.pressure)
        print("Iter CC Pressure: " + str(curr_cc_pressure))

        interval_massflow = self.calculate_injector_massflow(self.pressure, curr_cc_pressure)
        delta_system_mass = time_step * interval_massflow

        self.massflow = delta_system_mass

        if self.mass_liquid <= 0:
            self.mass_liquid = 0
            self.execute_vapor_phase_state()
            
            if self.mass_vapor <= 0:
                return False
            else:
                return True
            

        # Prior to accounting for nitrous vaporization to account 
        self.mass_liquid -= delta_system_mass 
        curr_system_mass = self.mass_liquid + self.mass_vapor

        true_liquid_mass = (self.volume  - (curr_system_mass / self.density_vapor))/ \
                (self.density_liquid**(-1) - self.density_vapor**(-1))

        mass_vaporized = self.mass_liquid - true_liquid_mass


        self.mass_vapor +=  mass_vaporized
        self.mass_liquid = true_liquid_mass
        
        self.execute_repressurization_cooling(mass_vaporized)

        if mass_vaporized < 0:
            return False # Negative vaporization mass means burnout


        return True


    def __init__(self, V_0, P_0, T_0, m_0, basic_ullage) -> None:
        
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





if __name__ == "__main__":
    exit_flag = False
    iter_count = 0

    NOS_model = NOS(0.04, 55, 298, 40, 0.15)

    sim_time = 0
    timestamps = []
    pressure_probe = []
    temp_probe = []
    mass_flow_probe = []
    TIME_STEP = 0.05


    while (iter_count < 10000 and not exit_flag):
        exit_flag = not NOS_model.execute_vapack(TIME_STEP)
        sim_time += TIME_STEP

        timestamps.append(sim_time)
        pressure_probe.append(NOS_model.pressure)
        temp_probe.append(NOS_model.temperature)
        mass_flow_probe.append(NOS_model.massflow)


        iter_count += 1

        if exit_flag:
            print("Terminated due to burnout")
            print("Simulation time: " + str(sim_time))

            if (DEBUG_VERBOSE):
                print(pressure_probe)
                print(temp_probe)
                print(mass_flow_probe)

            plt.plot(timestamps, pressure_probe, label = 'Pressure')
            plt.plot(timestamps, temp_probe, label = 'Temperature')
            plt.plot(timestamps, mass_flow_probe, label = 'Mass Flow')

            plt.legend()

            plt.show()

            
    if iter_count >= 998:
        print("Iteration Count exceeded")





        

