from time import perf_counter
from tank_blowdown_model import NOS_tank
import math
import matplotlib.pyplot as plt

import cantera as ct
import os 

PRINT_DEBUG_THERMOCHEM = True
PRINT_DEBUG_COMB_CHAMBER = True
PRINT_DEBUG_ENGINE = True

PRINT_DEBUG_SIM = True
PRINT_DEBUG_SIM_VERBOSE = True


BAR_TO_PA = 100000

R = 360

OUTPUT_FILE_PATH = "output.csv"

class N2O_HTPB_ThermochemistryModel:

    def __init__(self) -> None:
        self.gas_model = ct.Solution('Mevel2015-rocketry_modified.yaml')
    
    
    def sim_gas_mixture_combustion_temp(self, OF_ratio, temperature_K, pressure_Pa) -> float:

        # This requires a download of a NASA file and some alteration to
        # get it to compile correctly


        self.gas_model.TPY = temperature_K, pressure_Pa, f'N2O:{OF_ratio}, C4H6:1'

        self.gas_model.equilibrate('HP')

        return self.gas_model

class CombustionChamberModel:

    def __init__(self, D_0, L, _fuel_density) -> None: # All Dimensions in meters for now

        self.L = L
        self.D = D_0
        self.A = 0.25*math.pi*D_0**2
        self.cc_volume = self.A * self.L

        # regression constants
        self.a_regression_const = 0.397
        self.n_regression_const = 0.367
        self.m_regression_const = 1

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
        self.regression_rate = self.a_regression_const * (total_mass_flux ** self.n_regression_const) * (self.L ** self.m_regression_const) * 0.001 
        self.volumetric_regression_rate = self.regression_rate * (math.pi * self.D) * self.L # Perimeter calculation

        # Convert to massflow
        fuel_massflow = self.volumetric_regression_rate * self.fuel_density

        return fuel_massflow

    def sim_comubstion(self, oxidizer_massflow, delta_time):

        conv_flag = False
        self.fuel_mass_flux = 0 # Initial guess
        oxidizer_mass_flux = oxidizer_massflow / self.A

        if PRINT_DEBUG_COMB_CHAMBER:
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

            if print_flag and PRINT_DEBUG_COMB_CHAMBER:
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
        if PRINT_DEBUG_COMB_CHAMBER:
            print('Iteration regression rate: ' + str(self.regression_rate))
        self.A = 0.25*math.pi*self.D**2
        self.cc_volume = self.A * self.L
        self.fuel_mass_flux = self.fuel_mass_flux
        self.fuel_massflow  = self.fuel_mass_flux * self.A

        pass


class EngineModel:

    def __init__(self, params, area_ratio, throat_diam) -> None:
        self.tank_model = NOS_tank(*params)

        # density in kg/m^3, lengths in meters
        self.comb_chamber_model =\
                CombustionChamberModel(D_0=0.05, L=0.6, _fuel_density=1100)
        
        self.area_ratio = area_ratio
        self.combusted_gas = 'not-simulated'
        self.thermo_model = N2O_HTPB_ThermochemistryModel()
        self.cc_pressure = 101325 # Pa; initial pressure assumed to be atmospheric
        self.throat_diam = throat_diam
        self.A_throat = 0.25*math.pi*throat_diam**2

        self.cc_gas_mass = 0
        self.dP = 0
        self.dm = 0 

        self.pressure_build_time = 'not-evaluated'
        self.pressure_build_time_evaluated = False

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
        vapak_start = perf_counter()
        tank_status = self.tank_model.execute_vapack(delta_time, False, given_cc_pressure=self.cc_pressure/BAR_TO_PA)
        print('Vapak calculation time: ' + str(perf_counter() - vapak_start))

        if tank_status <= 0 or tank_status == 2:
            return -1 # Fuel exhausted


        # Calculate the mass flow rates into the combustion chamber
        ox_massflow = self.tank_model.massflow
        self.comb_chamber_model.sim_comubstion(oxidizer_massflow=ox_massflow, 
                delta_time=delta_time)
        fuel_massflow = self.comb_chamber_model.get_fuel_massflow()

        self.OF_ratio = ox_massflow/fuel_massflow
        print('OF Ratio: ' + str(self.OF_ratio))

        start_time = perf_counter()
        if self.combusted_gas == 'not-simulated':
            self.current_T = self.tank_model.temperature
        else:
            self.current_T = self.combusted_gas.T
        
        if (self.combusted_gas == 'not-simmed' or update_thermochem):
            self.combusted_gas = \
                    self.thermo_model.sim_gas_mixture_combustion_temp(\
                    OF_ratio=self.OF_ratio, temperature_K = 273.15, 
                    pressure_Pa = self.cc_pressure)
        if PRINT_DEBUG_SIM_VERBOSE:
            print('Cantera calculation time: ' + str(perf_counter() - start_time))

        print('Cantera combusted pressure: ' + str(self.combusted_gas.P))


        ### Calculate frozen flow nozzle parameters

        self.k = self.combusted_gas.cp / self.combusted_gas.cv
        print('Cantera calculated k: ' + str(self.k))

        # From eqn (3-22) RPE
        self.T_throat = 2 * self.combusted_gas.T / (self.k + 1)

        # See RPE chapter 3, pg 48 
        R = self.current_specific_R = 8314 / self.combusted_gas.mean_molecular_weight

        # RPE Equation (3-23); this is the mach 1 at the given temperature
        self.throat_velocity = math.sqrt(self.k*R*self.T_throat)

        # From flow mass continuity equation    
        self.throat_massflow = self.throat_velocity*self.A_throat*self.combusted_gas.density

        # trying out the adams bashford bullshit

        self.massflow_in = ox_massflow + fuel_massflow
        new_dm = ox_massflow + fuel_massflow - self.throat_massflow
        old_cc_mass = self.cc_gas_mass
        self.cc_gas_mass += delta_time * (new_dm * (3/2) - self.dm * (1/2))
        self.dm = new_dm


        print('CC Volume: ' + str(self.comb_chamber_model.cc_volume))

        new_dP = self.cc_pressure*((1/self.cc_gas_mass)*self.dm - 
                (1/self.comb_chamber_model.cc_volume)*self.comb_chamber_model.volumetric_regression_rate)

        print('dP: ' + str(new_dP))

        # trying out the adams bashford bullshit
        self.cc_pressure += delta_time * ((3/2)* new_dP - (1/2) * self.dP)
        self.dP = new_dP

        if self.dP < 0 and not self.pressure_build_time_evaluated:
            pass # the function has no internal concept of time

        print('CC Pressure: ' + str(self.cc_pressure))

        # From eqn (3-20) RPE
        # P_throat = self.combusted_gas.P * (2/(k + 1))**(k/(k - 1))
        P_throat = self.cc_pressure * (2/(self.k + 1))**(self.k/(self.k - 1))
        print('Throat pressure: ' + str(P_throat))

        # From a reverse-mogging of equation 3-25 of RPE
        P_exit = self.reverse_mog_exit_pressure(self.area_ratio, self.cc_pressure, self.k)
        
        print('Nozzle exit pressure: ' + str(P_exit))
        if PRINT_DEBUG_SIM_VERBOSE:
            print('Area ratio sanity check: ' + str(self.area_ratio_mog_equation((P_exit/self.cc_pressure), self.k)))


        # now just from isentropic equations:
        T_exit = self.T_throat * (P_exit/P_throat)**((self.k - 1)/self.k)

        
        self.velocity_exit = ((2*self.k/(self.k - 1))*R*self.combusted_gas.T * \
                        (1 - (P_exit/self.cc_pressure)**((self.k-1/self.k))))**(0.5)
        print('Nozzle exit velocity: ' + str(self.velocity_exit))

        # Now accounting for mass accumulation to calculate time-step dP: 
        # (will implement later)
        # CC_volume = self.comb_chamber_model.cc_volume
        # delta_P = combusted_gas.P * (1/(combusted_gas.density * CC_volume) * )

        thrust = (self.throat_massflow)*self.velocity_exit
        print('') # newline
        return thrust

if __name__ == '__main__':

    # Copy and pasted from the blowdown model code
    NOS_tank_params = [0.04, 55, 288, 40, 0.15]
    engine_model = EngineModel(NOS_tank_params, 4.8, 0.039385)
    # area_ratio = EngineModel.area_ratio_mog_equation(0.1, 1.2)
    # p_e = engine_model.reverse_mog_exit_pressure(area_ratio, 100000, 1.2)
    # print('remogged_ratio= ' + str(p_e/100000))
    


    sim_ok = True
    thrust_values = []
    cc_pressure_values = []
    timestamps = []
    massflow_out_values = []
    throat_temperature_values= []
    OF_values = []
    v_ex_values = []
    k_values = []
    massflow_in_values = []

    sim_time = 0
    SIM_STEP = 0.0125
    step_count = 0 

    # Setting this too low will throttle performance
    GRAPH_UPDATE_INTERVAL = 1.5

    fig1 = plt.figure(animated=True)

    ax1 = fig1.add_subplot(2,4,1)
    ax2 = fig1.add_subplot(2,4,2)
    ax3 = fig1.add_subplot(2,4,3)
    ax4 = fig1.add_subplot(2,4,4)
    ax5 = fig1.add_subplot(2,4,5)
    ax6 = fig1.add_subplot(2,4,6)
    ax7 = fig1.add_subplot(2,4,7)
    ax8 = fig1.add_subplot(2,4,8)
    
    plt.ion()

    plot1, = ax1.plot(timestamps, thrust_values)
    plot2, = ax2.plot(timestamps, cc_pressure_values)
    plot3, = ax3.plot(timestamps, massflow_out_values)
    plot4, = ax4.plot(timestamps, throat_temperature_values)
    plot5, = ax5.plot(timestamps, OF_values)
    plot6, = ax6.plot(timestamps, v_ex_values)
    plot7, = ax7.plot(timestamps, k_values)
    plot8, = ax8.plot(timestamps, massflow_in_values)
    

    fig1.canvas.draw()

    ax1.set_title('Thrust')
    ax2.set_title('CC Pressure')
    ax3.set_title('Massflow out')
    ax4.set_title('Throat temperature')
    ax5.set_title('OF ratio')
    ax6.set_title('Exhaust Velocity')
    ax7.set_title('K value')
    ax8.set_title('Massflow In')


    plt.show(block=False)

    axbackground = fig1.canvas.copy_from_bbox(ax1.bbox)
    ax2background = fig1.canvas.copy_from_bbox(ax2.bbox)

    ax1.autoscale_view()
    ax2.autoscale_view()
    
    prev_plot = perf_counter()

    while sim_ok:

        sim_thermochem = True
        if not ((step_count % 1) == 0):
            sim_thermochem = False

        sim_start = perf_counter()
        thrust_value = engine_model.sim_burn(SIM_STEP, update_thermochem=sim_thermochem)
        cc_pressure_value = engine_model.cc_pressure
        print('Simulation step compute time: ' + str(perf_counter() - sim_start))

        sim_time += SIM_STEP
        step_count += 1

        
        if perf_counter() - prev_plot > GRAPH_UPDATE_INTERVAL: 
            # ax1.clear()
            # ax2.clear()

            # ax1.set_title('Thrust')
            # ax2.set_title('CC Pressure')

            # ax1.plot(timestamps, thrust_values, color='red')
            # ax2.plot(timestamps, cc_pressure_values, color='red')

            plot1.set_data(timestamps, thrust_values)
            plot2.set_data(timestamps, cc_pressure_values)
            plot3.set_data(timestamps, massflow_out_values)
            plot4.set_data(timestamps, throat_temperature_values)
            plot5.set_data(timestamps, OF_values)
            plot6.set_data(timestamps, v_ex_values)
            plot7.set_data(timestamps, k_values)
            plot8.set_data(timestamps, massflow_in_values)

            ax1.relim()            
            ax2.relim()
            ax3.relim()            
            ax4.relim()
            ax5.relim()
            ax6.relim()
            ax7.relim()
            ax8.relim()

            ax1.autoscale_view()
            ax2.autoscale_view()
            ax3.autoscale_view()
            ax4.autoscale_view()
            ax5.autoscale_view()
            ax6.autoscale_view()
            ax7.autoscale_view()
            ax8.autoscale_view()

            # ax1.draw_artist(plot)
            # ax2.draw_artist(plot2)

            
            fig1.canvas.draw()
            prev_plot = perf_counter()

        fig1.canvas.flush_events()

        if thrust_value == -1:
            sim_ok = False # burnout

            plot1.set_data(timestamps, thrust_values)
            plot2.set_data(timestamps, cc_pressure_values)
            plot3.set_data(timestamps, massflow_out_values)
            plot4.set_data(timestamps, throat_temperature_values)
            plot5.set_data(timestamps, OF_values)
            plot6.set_data(timestamps, v_ex_values)
            plot7.set_data(timestamps, k_values)
            plot8.set_data(timestamps, massflow_in_values)

            ax1.relim()            
            ax2.relim()
            ax3.relim()            
            ax4.relim()
            ax5.relim()
            ax6.relim()
            ax7.relim()
            ax8.relim()

            ax1.autoscale_view()
            ax2.autoscale_view()
            ax3.autoscale_view()
            ax4.autoscale_view()
            ax5.autoscale_view()
            ax6.autoscale_view()
            ax7.autoscale_view()
            ax8.autoscale_view()

            # ax1.draw_artist(plot)
            # ax2.draw_artist(plot2)

            
            fig1.canvas.draw()
            
            with open(OUTPUT_FILE_PATH, 'w') as file:
                for idx in range(len(timestamps)):
                    file.write(str(timestamps[idx]) + ', ')
                    file.write(str(thrust_values[idx]) + ', ')
                    file.write(str(cc_pressure_values[idx]) + ', ')
                    file.write(str(massflow_out_values[idx]) + ', ')
                    file.write(str(throat_temperature_values[idx]) + ', ')
                    file.write(str(OF_values[idx]) + ', ')
                    file.write(str(k_values[idx]) + ', ')
                    file.write(str(massflow_in_values[idx]) + ', ')
                    file.write('\n')
            
            plt.show(block=True)

        else: 
            timestamps.append(sim_time)
            thrust_values.append(thrust_value)
            cc_pressure_values.append(cc_pressure_value)
            massflow_out_values.append(engine_model.throat_massflow)
            throat_temperature_values.append(engine_model.T_throat)
            OF_values.append(engine_model.OF_ratio)
            v_ex_values.append(engine_model.velocity_exit)
            k_values.append(engine_model.k)
            massflow_in_values.append(engine_model.massflow_in)
        
        # print('\n\n\n\n\n\n\n\n')

    # plt.show()

        # os.system('cls')
        print
        


