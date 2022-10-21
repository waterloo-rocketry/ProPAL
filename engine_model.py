from time import perf_counter
from tank_blowdown_model import NOS_tank
import math
import matplotlib.pyplot as plt

import cantera as ct

class N2O_HTPB_ThermochemistryModel:
    
    @staticmethod 
    def sim_gas_mixture_combustion_temp(OF_ratio, temperature_K, pressure_Pa) -> float:

        # This requires a download of a NASA file and some alteration to
        # get it to compile correctly

        gas_model = ct.Solution('Mevel2015-rocketry_modified.yaml')

        gas_model.TPY = temperature_K, pressure_Pa, f'N2O:{OF_ratio}, C4H6:1'

        gas_model.equilibrate('HP')

        return gas_model

class CombustionChamberModel:

    def __init__(self, D_0, L, _fuel_density) -> None: # All Dimensions in meters for now

        self.L = L
        self.D = D_0
        self.A = 0.25*math.pi*D_0**2
        self.cc_volume = self.A * self.L

        # regression constants
        self.a = 0.397
        self.n = 0.367
        self.m = 1

        self.fuel_density = _fuel_density 
        self.regression_rate = -1

        self.fuel_mass_flux = -1
        self.fuel_massflow = -1

        pass
        

    def get_fuel_mass_flux(self):
        if self.fuel_mass_flux == -1:
            raise RuntimeError('Fuel mass flux has not been simulated yet')
        
        return self.fuel_mass_flux

    def get_fuel_massflow(self):
        if self.fuel_massflow == -1:
            raise RuntimeError('Fuel massflow has not been simulated yet')
        
        return self.fuel_massflow

    def sim_fuel_regression_massflow(self, total_mass_flux, delta_time):
        # Constants are tuned for mm/s
        self.regression_rate = self.a * (total_mass_flux ** self.n) * (self.L ** self.m) * 0.001 
        fuel_volume_regressed_during_iter = (self.regression_rate * delta_time) * (math.pi * self.D) * self.L

        # Convert back to rate
        fuel_massflow = fuel_volume_regressed_during_iter * self.fuel_density / delta_time

        return fuel_massflow

    def sim_comubstion(self, oxidizer_massflow, delta_time):

        conv_flag = False
        self.fuel_mass_flux = 0 # Initial guess
        oxidizer_mass_flux = oxidizer_massflow / self.A

        print('Oxidizer massflux = ' + str(oxidizer_mass_flux))
        print('Port area = ' + str(self.A))

        iter_count = 0

        while not conv_flag:
            print_flag = ((iter_count % 1) == 0) and False
            total_mass_flux = oxidizer_mass_flux + self.fuel_mass_flux

            if print_flag:
                print('Total mass flux: ' + str(total_mass_flux))
                print('')
            
            recalc_fuel_mass_flux = \
                    (self.sim_fuel_regression_massflow(total_mass_flux=total_mass_flux,
                    delta_time=delta_time) / self.A)

            conv_criteria = (abs(self.fuel_mass_flux - recalc_fuel_mass_flux))/(self.fuel_mass_flux + 1e-9)

            if print_flag:
                print('Recalculated flux: ' + str(recalc_fuel_mass_flux))
                print('Convergence criteria: ' + str(conv_criteria))

            if conv_criteria  < 0.001:
                conv_flag = True                
                print('Converged fuel flux: ' + str(recalc_fuel_mass_flux) + ' Converged fuel massflow: ' + str(recalc_fuel_mass_flux * self.A))
            
            self.fuel_mass_flux = recalc_fuel_mass_flux

            iter_count += 1
            if iter_count >= 10000:
                raise RuntimeError('Combustion Chamber fuel flux failed to converge')


        self.D += 2*self.regression_rate*delta_time
        print('Iteration regression rate: ' + str(self.regression_rate))
        self.A = 0.25*math.pi*self.D**2
        self.cc_volume = self.A * self.L
        self.fuel_mass_flux = self.fuel_mass_flux
        self.fuel_massflow  = self.fuel_mass_flux * self.A

        pass


class EngineModel:

    def __init__(self, params, area_ratio) -> None:
        self.tank_model = NOS_tank(*params)

        # density in kg/m^3, lengths in meters
        self.comb_chamber_model =\
                CombustionChamberModel(D_0=0.05, L=0.6, _fuel_density=1100)
        
        self.area_ratio = area_ratio
        self.combusted_gas = 'not-simulated'

    @classmethod
    def reverse_mog_exit_pressure(cls, area_ratio, p1, k):
        conv_flag = False

        # p_exit_guess = p1 - 1e-5  # Pa
        # p_increment = -100

        p_exit_guess = 0 + 1e-5
        p_increment = 1000

        iter_count = 0

        while not conv_flag:
            p_exit_guess += p_increment

            
            # print('Pressure guess: ' + str(p_exit_guess))
            # print('Pressure ratio: ' + str(p_exit_guess/p1))
            resulting_area_ratio = cls.area_ratio_mog_equation(p_exit_guess/p1, k)
            # print('Computed area ratio: ' + str(resulting_area_ratio))

            if p_increment > 0 and resulting_area_ratio < area_ratio:
                p_increment = -math.sqrt(p_increment)
            elif p_increment < 0 and resulting_area_ratio > area_ratio:
                p_increment = math.sqrt(abs(p_increment))

            iter_count += 1

            if (abs(resulting_area_ratio - area_ratio))/area_ratio < 0.0001:
                # print('Exit pressure converged in ' + str(iter_count) + ' cycles')
                conv_flag = True
                # print(p_increment)

            if iter_count >= 10000:
                raise RuntimeError('Mogged exit pressure failed to converge')



        return p_exit_guess


    @staticmethod
    def area_ratio_mog_equation(p_ratio, k):
        return (((k + 1.0)/2.0)**(1.0/(k - 1.0)) * \
                (p_ratio)**(1.0/k) * \
                ( ((k + 1.0)/(k - 1.0)) * (1.0 - (p_ratio)**((k-1.0)/(k))) )**(0.5))**(-1.0)

    @classmethod
    def test_area_ratio_mog_equation(cls):
        x = []
        mog = []
        k = 1.2
        for i in range (10, 100):
            x.append(i)
            mog.append(cls.area_ratio_mog_equation(i**(-1.0), k))
        
        plt.plot(x, mog)
        plt.show()



    def sim_burn(self, delta_time, update_thermochem = True):
        tank_status = self.tank_model.execute_vapack(delta_time)
        if tank_status <= 0 or tank_status == 2:
            return -1 # Fuel exhausted


        # Calculate the mass flow rates into the combustion chamber
        ox_massflow = self.tank_model.massflow
        self.comb_chamber_model.sim_comubstion(oxidizer_massflow=ox_massflow, 
                delta_time=delta_time)
        fuel_massflow = self.comb_chamber_model.get_fuel_massflow()

        OF_ratio = ox_massflow/fuel_massflow
        print('OF Ratio: ' + str(OF_ratio))

        bar_to_Pa = 100000

        start_time = perf_counter()
        if (self.combusted_gas == 'not-simmed' or update_thermochem):
            self.combusted_gas = \
                    N2O_HTPB_ThermochemistryModel.sim_gas_mixture_combustion_temp(\
                    OF_ratio=OF_ratio, temperature_K = 298, 
                    pressure_Pa = bar_to_Pa*(self.tank_model.pressure/2))
        print('Cantera calculation time: ' + str(perf_counter() - start_time))
        print('Cantera combusted pressure: ' + str(self.combusted_gas.P))


        ### Calculate frozen flow nozzle parameters

        k = self.combusted_gas.cp / self.combusted_gas.cv
        print('Cantera calculated k: ' + str(k))

        # From eqn (3-22) RPE
        T_throat = 2 * self.combusted_gas.T / (k + 1)

        # From eqn (3-20) RPE
        P_throat = self.combusted_gas.P * (2/(k + 1))**(k/(k - 1))
        print('Throat pressure: ' + str(P_throat))

        # From a reverse-mogging of equation 3-25 o RPE
        P_exit = self.reverse_mog_exit_pressure(self.area_ratio, self.combusted_gas.P, k)
        print('Nozzle exit pressure: ' + str(P_exit))
        print('Area ratio sanity check: ' + str(self.area_ratio_mog_equation((P_exit/self.combusted_gas.P), k)))


        # now just from isentropic equations:
        T_exit = T_throat * (P_exit/P_throat)**((k - 1)/k)

        R = 360
        velocity_exit = ((2*k/(k - 1))*R*self.combusted_gas.T * \
                        (1 - (P_exit/self.combusted_gas.P)**((k-1/k))))**(0.5)
        print('Nozzle exit velocity: ' + str(velocity_exit))

        # Now accounting for mass accumulation to calculate time-step dP: 
        # (will implement later)
        # CC_volume = self.comb_chamber_model.cc_volume
        # delta_P = combusted_gas.P * (1/(combusted_gas.density * CC_volume) * )

        thrust = (ox_massflow + fuel_massflow)*velocity_exit
        print('') # newline
        return thrust

if __name__ == '__main__':

    # Copy and pasted from the blowdown model code
    NOS_tank_params = [0.04, 55, 288, 40, 0.15]
    engine_model = EngineModel(NOS_tank_params, 4.8)
    # area_ratio = EngineModel.area_ratio_mog_equation(0.1, 1.2)
    # p_e = engine_model.reverse_mog_exit_pressure(area_ratio, 100000, 1.2)
    # print('remogged_ratio= ' + str(p_e/100000))
    


    sim_ok = True
    thrust_values = []
    timestamps = []

    sim_time = 0
    SIM_STEP = 0.025
    step_count = 0 

    while sim_ok:

        sim_thermochem = True
        if not ((step_count % 1) == 0):
            sim_thermochem = False

        thrust_value = engine_model.sim_burn(SIM_STEP, update_thermochem=sim_thermochem)
        sim_time += SIM_STEP
        step_count += 1

        if thrust_value == -1:
            sim_ok = False # burnout
        else:
            thrust_values.append(thrust_value)
            timestamps.append(sim_time)

    plt.plot(timestamps, thrust_values)
    plt.title('Thrust curve')
    plt.show()

