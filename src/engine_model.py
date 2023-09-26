_imports = 0 # Outline checkpoint

from copy import deepcopy
from time import perf_counter
from ruamel.yaml import YAML

import math
import matplotlib.pyplot as plt
import os 

from tank_blowdown_model import NOS_tank
from OptimizationTemp.simulation_generator import LookUpTable
from ThermochemWrappers import N2O_HTPB_ThermochemistryModel

PRINT_DEBUG_THERMOCHEM = False
PRINT_DEBUG_COMB_CHAMBER = False
PRINT_DEBUG_ENGINE = False
PRINT_DEBUG_SIM = False
PRINT_DEBUG_SIM_VERBOSE = False


BAR_TO_PA = 100000

R = 360

OUTPUT_FILE_PATH = "output.csv"

try:
    with open("input.yaml", "r") as file:
        yaml = YAML(typ="safe")
        config = yaml.load(file)
except FileNotFoundError:
    print("input.yaml file required")
    raise SystemExit

PRINT_DEBUG_THERMOCHEM |= config["flags"]["debug"]
PRINT_DEBUG_COMB_CHAMBER |= config["flags"]["debug"]
PRINT_DEBUG_ENGINE |= config["flags"]["debug"]
PRINT_DEBUG_SIM |= config["flags"]["debug"]

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

        # For third order Adams-Bashforth integration
        self.previous_regression_rate = 0
        self.regression_rate = 0
        self.new_regression_rate = 0

        self.previous_mass_flux = -1
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

    def sim_fuel_regression_massflow(self, total_mass_flux):
        # Constants are tuned for mm/s
        self.new_regression_rate = self.a_regression_const * (total_mass_flux ** self.n_regression_const) * (self.L ** self.m_regression_const) * 0.001 
        self.volumetric_regression_rate = self.new_regression_rate * (math.pi * self.D) * self.L # Perimeter calculation

        # Convert to massflow
        fuel_massflow = self.volumetric_regression_rate * self.fuel_density

        return fuel_massflow

    def sim_comubstion(self, oxidizer_massflow, delta_time):

        conv_flag = False
        self.fuel_mass_flux = 1e-9 # Initial guess
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
                    (self.sim_fuel_regression_massflow(total_mass_flux=total_mass_flux) / self.A)

            conv_criteria = (abs(self.fuel_mass_flux - recalc_fuel_mass_flux))/(self.fuel_mass_flux)

            if print_flag and PRINT_DEBUG_COMB_CHAMBER:
                print('Recalculated flux: ' + str(recalc_fuel_mass_flux))
                print('Convergence criteria: ' + str(conv_criteria))

            if conv_criteria  < 0.0001:
                conv_flag = True     
                if PRINT_DEBUG_COMB_CHAMBER:           
                    print('Converged fuel flux: ' + str(recalc_fuel_mass_flux) +\
                         ' Converged fuel massflow: ' + str(recalc_fuel_mass_flux * self.A))
            
            self.fuel_mass_flux = recalc_fuel_mass_flux

            iter_count += 1
            if iter_count >= 100000:
                raise RuntimeError('Combustion Chamber fuel flux failed to converge')


        self.D += 2 * delta_time * ((23/12)*self.new_regression_rate - \
                (16/12) * self.regression_rate + (5/12)*self.previous_regression_rate)

        self.previous_regression_rate = self.regression_rate
        self.regression_rate = self.new_regression_rate
        if PRINT_DEBUG_COMB_CHAMBER:
            print('Iteration regression rate: ' + str(self.regression_rate))
        self.A = 0.25*math.pi*self.D**2
        self.cc_volume = self.A * self.L
        self.fuel_massflow  = self.fuel_mass_flux * self.A

        pass


class EngineModel:

    def __init__(self, nos_tank_params, area_ratio, throat_diam, combustion_efficiency, use_external_tc_model=False) -> None:

        if type(nos_tank_params) == NOS_tank:
            self.tank_model = nos_tank_params
        else:
            self.tank_model = NOS_tank(*nos_tank_params)

        # density in kg/m^3, lengths in meters
        self.comb_chamber_model =\
                CombustionChamberModel(D_0=0.05, L=0.6, _fuel_density=1100)
        
        self.area_ratio = area_ratio
        self.combusted_gas = 'not-simulated'

        if not use_external_tc_model:
            self.thermo_model = N2O_HTPB_ThermochemistryModel()
        else:
            self.thermo_model = None


        self.cc_pressure = 101325 # Pa; initial pressure assumed to be atmospheric
        self.throat_diam = throat_diam
        self.A_throat = 0.25*math.pi*throat_diam**2

        self.combustion_efficiency = combustion_efficiency

        self.cc_gas_mass = 0

        self.previous_dP = 0
        self.dP = 0

        self.previous_dm = 0
        self.dm = 0 


        self.elapsed_time = 0
        self.pressure_build_time = 'not-evaluated'
        self.pressure_build_time_evaluated = False
        self.thurst_at_pressure_peak = 'not-evaluated'

        self.choked = False

        self.iter_count = 0


    @classmethod
    def solve_exit_pressure(cls, area_rat00io, p1, k, choked = True):
        '''
        Doesn't work yet. 
        '''
        from scipy.optimize import fsolve

        fsolve()



    @classmethod
    def reverse_mog_exit_pressure(cls, area_ratio, p1, k, choked = True, pressure_ratio_only = False):
        '''
        This hand-written numerical solver was written before I knew what I was doing. It exists for 
        reference, but a proper solver from a solution library will be used as soon as it can be
        implemented successfully. 
        '''

        conv_flag = False

        # p_exit_guess = p1 - 1e-5  # Pa
        # p_increment = -100

        if choked:
            p_exit_guess = 0 + 1e-6
            p_increment = 10000
        else: # if flow is not choked, it will try to find the higher (subsonic) solution
            p_exit_guess = p1 - 1e-6
            p_increment = -10000

        iter_count = 0

        while not conv_flag:
            

            
            # print('Pressure guess: ' + str(p_exit_guess))
            # print('Pressure ratio: ' + str(p_exit_guess/p1))
            resulting_area_ratio = cls.area_ratio_mog_equation(p_exit_guess/p1, k)
            # print('Computed area ratio: ' + str(resulting_area_ratio))
            # print(area_ratio)



            iter_count += 1

            if (abs(resulting_area_ratio - area_ratio))/area_ratio < 0.00025:
                # print('Exit pressure converged in ' + str(iter_count) + ' cycles')
                conv_flag = True
                # print(p_increment)

            if iter_count >= 100000:
                raise RuntimeError('Mogged exit pressure failed to converge')

            
            if choked:
                if p_increment > 0 and resulting_area_ratio < area_ratio:
                    p_increment = -math.sqrt(p_increment)
                elif p_increment < 0 and resulting_area_ratio > area_ratio:
                    p_increment = math.sqrt(abs(p_increment))
            else:
                if p_increment > 0 and resulting_area_ratio > area_ratio:
                    p_increment = -math.sqrt(p_increment)
                elif p_increment < 0 and resulting_area_ratio < area_ratio:
                    p_increment = math.sqrt(abs(p_increment))


            p_exit_guess += p_increment

           


        if pressure_ratio_only:
            return p_exit_guess/p1
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



    def sim_burn(self, delta_time, update_thermochem = True, external_tc_model = None, suppress_prints = False):


        vapak_start = perf_counter()
        tank_status = self.tank_model.execute_vapack(delta_time, suppress_prints= True, given_cc_pressure=self.cc_pressure/BAR_TO_PA)
        if PRINT_DEBUG_SIM_VERBOSE:
            print('Vapak calculation time: ' + str(perf_counter() - vapak_start))

        if tank_status <= 0 or tank_status == 2:
            return -1 # Fuel exhausted


        # Calculate the mass flow rates into the combustion chamber
        ox_massflow = self.tank_model.massflow
        self.comb_chamber_model.sim_comubstion(oxidizer_massflow=ox_massflow, 
                delta_time=delta_time)
        fuel_massflow = self.comb_chamber_model.get_fuel_massflow()

        self.OF_ratio = ox_massflow/fuel_massflow

        if PRINT_DEBUG_ENGINE:
            print('OF Ratio: ' + str(self.OF_ratio))

        start_time = perf_counter()
        

        if self.thermo_model:
            sim_tc_model = self.thermo_model
        else:
            sim_tc_model = external_tc_model

        if (update_thermochem):
            combusted_gas = \
                    sim_tc_model.sim_gas_mixture_combustion_temp(\
                    OF_ratio=self.OF_ratio, temperature_K = 273.31, # Current_T not used yet
                    pressure_Pa = self.cc_pressure)
        if PRINT_DEBUG_SIM_VERBOSE:
            print('Cantera calculation time: ' + str(perf_counter() - start_time))

        if PRINT_DEBUG_THERMOCHEM:
            print('Cantera combusted pressure: ' + str(combusted_gas.P))


        ### Calculate frozen flow nozzle parameters

        self.k = combusted_gas.cp / combusted_gas.cv
        k = self.k

        if PRINT_DEBUG_ENGINE:
            print('Cantera calculated k: ' + str(self.k))



        # See RPE chapter 3, pg 48 
        R = self.current_specific_R = 8314 / combusted_gas.mean_molecular_weight

        # It needs to be determined if the flow is even sonic:
        critical_upstream_pressure = 101325 * (2/(k + 1))**(-k/(k - 1)) # From wikipedia, choked flow condition. 
        if self.cc_pressure < critical_upstream_pressure:
            # https://www.engineersedge.com/pressure,045vessel/gas_discharge_rate_14170.htm
            # This blessed link gives us the salvation we need
            DISCHARGE_COEFF = 0.72
            P_a = 101325
            P_ratio = P_a/self.cc_pressure
            self.throat_massflow = DISCHARGE_COEFF * self.A_throat * \
                    math.sqrt( 2*self.cc_pressure*combusted_gas.density*(k/(k - 1)))* \
                    (P_ratio**(2/k) - P_ratio**((k+1)/k))
            self.choked = False
        else:
            # From flow mass continuity equation    
            # RPE Equation (3-23); this is the mach 1 at the given temperature
            # self.throat_velocity = math.sqrt(self.k*R*self.T_throat)
            # self.throat_massflow = self.throat_velocity*self.A_throat*self.combusted_gas.density

            # From RPE eqn 3-24
            self.throat_massflow = self.A_throat * self.cc_pressure * k * \
                    math.sqrt(((2/(k + 1))**((k + 1)/(k - 1))))/(math.sqrt(k * R * combusted_gas.T))

            self.choked = True
            

        # From eqn (3-22) RPE
        self.T_throat = 2 * combusted_gas.T / (self.k + 1)




        # This is assuming an always choked flow: not necessarily true


        # trying out the adams bashforth multistep integration
        self.massflow_in = ox_massflow + fuel_massflow

        new_dm = ox_massflow + fuel_massflow - self.throat_massflow
        # differentially_mogged_dm = self.cc_gas_mass*((1/self.cc_pressure)*self.dP + 
        #         (1/self.comb_chamber_model.cc_volume)*self.comb_chamber_model.volumetric_regression_rate)
        # new_dm = differentially_mogged_dm
        # # back-calculated
        # self.throat_massflow = self.massflow_in - new_dm

        # # trying out the adams bashforth multistep integration (3rd order)

        delta_m_second_order = delta_time * ((3/2) * new_dm - (1/2) * self.dm)
        delta_m_third_order = delta_time * ((23/12) * new_dm - (16/12) * self.dm + (5/12)*self.previous_dm)

        final_step_delta_m = 3*delta_m_third_order - 2*delta_m_second_order

        # Use higher-order estimates as they become available due to previous calculations
        if self.iter_count == 0:
            final_step_delta_m = delta_time * new_dm
        elif self.iter_count == 1:
            final_step_delta_m = delta_m_second_order
        else:
            final_step_delta_m = delta_m_third_order

        # Trying some richardson extrapolation bullshit, idk if this will even work
        self.cc_gas_mass += delta_time * new_dm

        self.numerical_error = abs(delta_m_third_order - delta_m_second_order)

        self.previous_dm = self.dm
        self.dm = new_dm

        if PRINT_DEBUG_COMB_CHAMBER:
            print('CC Volume: ' + str(self.comb_chamber_model.cc_volume))

        new_dP = self.cc_pressure*((1/self.cc_gas_mass)*self.dm - 
                (1/self.comb_chamber_model.cc_volume)*self.comb_chamber_model.volumetric_regression_rate)

        if PRINT_DEBUG_ENGINE:
            print('dP: ' + str(new_dP))

        # Trying out the adams bashforth multistep integration (3rd order)

        # Trying some richardson extrapolation bullshit, idk if this will even work
        delta_p_second_order = delta_time * ((3/2) * new_dP - (1/2) * self.dP)
        delta_p_third_order = delta_time * ((23/12) * new_dP - (16/12) * self.dP + (5/12)*self.previous_dP)

        final_step_delta_P = 3*delta_p_third_order - 2*delta_p_second_order

        # Use higher-order estimates as they become available due to previous calculations
        if self.iter_count == 0:
            final_step_delta_P = delta_time * new_dP
        elif self.iter_count == 1:
            final_step_delta_P = delta_p_second_order
        else:
            final_step_delta_P = delta_p_third_order

        self.cc_pressure += delta_time * new_dP
        self.previous_dP = self.dP
        self.dP = new_dP

        if PRINT_DEBUG_ENGINE:
            print('CC Pressure: ' + str(self.cc_pressure))

        # From eqn (3-20) RPE
        # P_throat = self.combusted_gas.P * (2/(k + 1))**(k/(k - 1))
        P_throat = self.cc_pressure * (2/(self.k + 1))**(self.k/(self.k - 1))
        if PRINT_DEBUG_ENGINE:
            print('Throat pressure: ' + str(P_throat))

        # From a reverse-mogging of equation 3-25 of RPE
        P_exit = self.reverse_mog_exit_pressure(self.area_ratio, self.cc_pressure, self.k, choked=self.choked)
        if PRINT_DEBUG_ENGINE:
            print('Nozzle exit pressure: ' + str(P_exit))
        if PRINT_DEBUG_SIM_VERBOSE:
            print('Area ratio sanity check: ' + str(self.area_ratio_mog_equation((P_exit/self.cc_pressure), self.k)))

        
        self.velocity_exit = ((2*self.k/(self.k - 1))*R*combusted_gas.T * \
                        (1 - (P_exit/self.cc_pressure)**((self.k-1/self.k))))**(0.5)
        if PRINT_DEBUG_ENGINE:
            print('Nozzle exit velocity: ' + str(self.velocity_exit))

        self.thrust = (self.throat_massflow)*self.velocity_exit * self.combustion_efficiency

        if self.dP < 0 and not self.pressure_build_time_evaluated:
            self.pressure_build_time = self.elapsed_time
            self.pressure_build_time_evaluated = True
            self.thurst_at_pressure_peak = self.thrust

        print('') # newline

        # For tracking how numerical methods are implemented
        self.elapsed_time += delta_time
        self.iter_count += 1

        return self.thrust

class MockLookupThermoModel():
    
    def __init__(self) -> None:
        self.lookup_table = \
                LookUpTable(os.path.join(os.getcwd(), "src", "OptimizationTemp", "SG-TargetedTableOF_ratio" +\
                "-temperature_K-pressure_Pa---0.1-15-0.25-273.3-273.3-0.1-1-10000000.0-25000.csv"))

    def sim_gas_mixture_combustion_temp(self, OF_ratio, temperature_K, pressure_Pa) -> float:
        return self.lookup_table.lookup(OF_ratio, temperature_K, pressure_Pa)



class HybridBurnSimulator:
    @staticmethod
    def sim_full_burn():

        updated_tank_params = [0.0235, -99, 302.4, -99, 0.1]

        # Copy and pasted from the blowdown model code

        # V_0, P_0, T_0, m_0, basic_ullage
        # NOS_tank_params = [0.04, 55, 288, 40, 0.15]
        NOS_tank_params = [config["oxidiser"]["V_0"], config["oxidiser"]["P_0"], config["oxidiser"]["T_0"], \
        config["oxidiser"]["m_0"], config["oxidiser"]["ullage"]]

        model_creation_start = perf_counter()

        sim_tank_model = NOS_tank(*updated_tank_params)
        
        engine_model = EngineModel(sim_tank_model, config["engine"]["area_ratio"], config["engine"]["throat_d"], config["engine"]["combust_e"], use_external_tc_model=config["flags"]["external_tc"])

        model_creation_end = perf_counter()

        print("Engine model creation time: " + str(- model_creation_start + model_creation_end))

        USE_REAL_THERMO_MODEL = False

        if USE_REAL_THERMO_MODEL:
            thermo_model = N2O_HTPB_ThermochemistryModel()
        else:
            thermo_model = MockLookupThermoModel()

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

        numerical_error_values = []
        step_error_values = []
        step_time_values = []

        sim_time = 0
        INITIAL_SIM_STEP = 5e-5
        NUMERICAL_EPSILON = 1e-6
        step_count = 0 

        DEBUG_STEP_LIMIT = -99 # For halting the program a set amount of iterations in (negative values will run the simulation normally)

        # Setting this too low will throttle performance
        GRAPH_UPDATE_INTERVAL = 2.5



        fig1 = plt.figure(animated=True, figsize = (16,8))
        plt.subplots_adjust(wspace=0.4, hspace=0.4, bottom=0.05, left=0.05, top=0.95, right=0.95) 

        ax1 = fig1.add_subplot(3,5,1)
        ax2 = fig1.add_subplot(3,5,2)
        ax3 = fig1.add_subplot(3,5,3)
        ax4 = fig1.add_subplot(3,5,4)
        ax5 = fig1.add_subplot(3,5,5)
        ax6 = fig1.add_subplot(3,5,6)
        ax7 = fig1.add_subplot(3,5,7)
        ax8 = fig1.add_subplot(3,5,8)
        ax9 = fig1.add_subplot(3,5,9)

        ax10 = fig1.add_subplot(3,5,10) # for print-outs

        ax11 = fig1.add_subplot(3,5,11)
        ax12 = fig1.add_subplot(3,5,12)

        ax13 = fig1.add_subplot(3,5,13)
        ax14 = fig1.add_subplot(3,5,14)


        
        plt.ion()
        

        plot1, = ax1.plot(timestamps, thrust_values)
        ax1.text(0.1,0.9, f'Efficiency: {engine_model.combustion_efficiency:.2f}', 
                            ha = 'left', transform=ax1.transAxes, fontsize = 'xx-small')

        plot2, = ax2.plot(timestamps, cc_pressure_values)
        plot3, = ax3.plot(timestamps, massflow_out_values)
        plot4, = ax4.plot(timestamps, throat_temperature_values)
        plot5, = ax5.plot(timestamps, OF_values)
        plot6, = ax6.plot(timestamps, v_ex_values)
        plot7, = ax7.plot(timestamps, k_values)
        plot8, = ax8.plot(timestamps, massflow_in_values)
        plot9, = ax9.plot(timestamps, numerical_error_values)

        plot11, = ax11.plot(timestamps, step_error_values)
        plot12, = ax12.plot(timestamps, step_time_values)
        

        fig1.canvas.draw()

        ax1.set_title('Thrust')
        ax2.set_title('CC Pressure')
        ax3.set_title('Massflow out')
        ax4.set_title('Throat temperature')
        ax5.set_title('OF ratio')
        ax6.set_title('Exhaust Velocity')
        ax7.set_title('K value')
        ax8.set_title('Massflow In')

        ax9.set_title('Numerical Error')

        ax11.set_title('Step Error')
        ax12.set_title('Time Step')

        axs = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax10,ax9, ax11,ax12]

        for ax in axs:
            plt.sca(ax)

            plt.xticks(fontsize = 'small')
            plt.yticks(fontsize = 'small')
            # ax.yaxis.set_major_formatter(FormatStrFormatter('%2.1e'))
        
        plt.sca(ax10)
        plt.tick_params(
            axis='both',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom=False,      # ticks along the bottom edge are off
            left = False,
            labelbottom=False,
            labelleft=False
            )
    


        plt.show(block=False)
        
        prev_plot = perf_counter()

        ADAPTIVE_STEP = False

        # Absolute and relative errors for the variable being measured (currently pressure in Pa)
        # The error scale is thus atol + P*rtol
        STEP_RTOL = 1e-6
        STEP_ATOL = 1

        step_size = INITIAL_SIM_STEP


        while sim_ok:

            sim_thermochem = True
            if not ((step_count % 1) == 0):
                sim_thermochem = False

            

            step_start_time = perf_counter()
            step_error = 1e12

            if ADAPTIVE_STEP:
                
                iter_abort_flag = False
                adaptive_step_iters = 0
                previous_step_size = step_size

                while step_error > NUMERICAL_EPSILON and not iter_abort_flag:

                    full_step_model = deepcopy(engine_model)
                    half_step_model = deepcopy(engine_model)

                    full_step_model.sim_burn(step_size, update_thermochem=sim_thermochem,
                             external_tc_model=thermo_model)

                    half_step_model.sim_burn(step_size/2, update_thermochem=sim_thermochem, 
                            external_tc_model=thermo_model)
                    half_step_model.sim_burn(step_size/2, update_thermochem=sim_thermochem, 
                            external_tc_model=thermo_model)
                    
                    full_step_pressure = full_step_model.cc_pressure
                    half_step_pressure = half_step_model.cc_pressure

                    # Most of this numerical step size code follow the approach suggested by 
                    # Press et al (2007) in Numerical Recipes 3e, pg 913-915

                    truncation_error = full_step_pressure - half_step_pressure

                    err_scale = STEP_ATOL + STEP_RTOL*full_step_pressure

                    step_error = abs(truncation_error/err_scale)
                    if step_error == 0:
                        break

                    TARGET_ERROR = 1 # err_scale/err_scale
                    ITER_SAFETY_FACTOR = 0.98

                    # Taken as a modification of Numerical Recipes eqn (17.2.10)
                    # [The original approach was for Runge Kutta 5th order, the existing 
                    # solution is 3rd order, hence the use of 1/3 instead of 1/5 in the exponent
                    # for more details, read the textbook]
                    
                    # This was made up by me, due to the overall dubious nature of this
                    # code; in order to prevent oscillation. This probably slows down the code
                    # significantly
                    COARSENING_SAFTY_FACTOR = 0.975

                    if step_error < COARSENING_SAFTY_FACTOR * TARGET_ERROR: 
                        # Step size can be coarsened, and then the function can exit:
                        step_size = ITER_SAFETY_FACTOR * step_size *\
                                (TARGET_ERROR/step_error)**(1/3)
                        iter_abort_flag = True
                        print(adaptive_step_iters)
                    elif step_error > TARGET_ERROR:
                        # Error too great, the step size must be coarsened
                        step_size = ITER_SAFETY_FACTOR * step_size *\
                                (TARGET_ERROR/step_error)**(1/3)
                    else:
                        # No modification of step size required. 
                        iter_abort_flag = True
                        print(adaptive_step_iters)



                    if adaptive_step_iters >= 10:
                        print('WARNING: Iteration count (10) for adaptive step exceeded!')
                        iter_abort_flag = True
                    
                    if step_size >= previous_step_size * 10 or step_size <= previous_step_size * 5:
                        iter_abort_flag = True

                    adaptive_step_iters += 1
                    # iter_abort_flag = True

                engine_model = deepcopy(half_step_model)

                thrust_value = engine_model.thrust
                cc_pressure_value = engine_model.cc_pressure

            else:

                thrust_value = engine_model.sim_burn(INITIAL_SIM_STEP, update_thermochem=sim_thermochem, external_tc_model=thermo_model)
                cc_pressure_value = engine_model.cc_pressure

            if PRINT_DEBUG_SIM_VERBOSE:
                print('Simulation step compute time: ' + str(perf_counter() - step_start_time))

            sim_time += step_size
            step_count += 1

            
            if perf_counter() - prev_plot > GRAPH_UPDATE_INTERVAL: 
                plot1.set_data(timestamps, thrust_values)
                plot2.set_data(timestamps, cc_pressure_values)
                plot3.set_data(timestamps, massflow_out_values)
                plot4.set_data(timestamps, throat_temperature_values)
                plot5.set_data(timestamps, OF_values)
                plot6.set_data(timestamps, v_ex_values)
                plot7.set_data(timestamps, k_values)
                plot8.set_data(timestamps, massflow_in_values)
                plot9.set_data(timestamps, numerical_error_values)

                plot11.set_data(timestamps, step_error_values)
                plot12.set_data(timestamps, step_time_values)
                

                for ax in axs:
                    ax.relim()
                    ax.autoscale_view()

                fig1.canvas.draw()
                
                with open(OUTPUT_FILE_PATH, 'w') as file:
                    file.write('Time (s), Thrust (N), CC Pressure (Pa), '+\
                        'mdot_out (kg/s), T_t (K), OF (dimless), k (dimless), mdot_in (kg/s)')

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

                ax10.clear()
                ax10.text(0.9,0.7, f'CC Pres: {engine_model.cc_pressure:.3f}', 
                            ha = 'right', transform=ax10.transAxes, fontsize = 'small')
                ax10.text(0.9,0.6, f'Choked?: {engine_model.choked}', 
                            ha = 'right', transform=ax10.transAxes, fontsize = 'small')

                if engine_model.pressure_build_time_evaluated:
                    ax10.text(0.9,0.9, f'Pres rise time: {engine_model.pressure_build_time:.3f}', 
                            ha = 'right', transform=ax10.transAxes, fontsize = 'small')
                    ax10.text(0.9,0.8, f'Peak thrust: {engine_model.thurst_at_pressure_peak:.3f}', 
                            ha = 'right', transform=ax10.transAxes, fontsize = 'small')
                    text_display_flag = True

                
                fig1.canvas.draw()
                prev_plot = perf_counter()

            fig1.canvas.flush_events()

            if thrust_value == -1 or step_count == DEBUG_STEP_LIMIT: # Stop the program to analyze debug
                sim_ok = False # burnout

                plot1.set_data(timestamps, thrust_values)
                plot2.set_data(timestamps, cc_pressure_values)
                plot3.set_data(timestamps, massflow_out_values)
                plot4.set_data(timestamps, throat_temperature_values)
                plot5.set_data(timestamps, OF_values)
                plot6.set_data(timestamps, v_ex_values)
                plot7.set_data(timestamps, k_values)
                plot8.set_data(timestamps, massflow_in_values)
                plot9.set_data(timestamps, numerical_error_values)

                plot11.set_data(timestamps, step_error_values)
                plot12.set_data(timestamps, step_time_values)
                

                for ax in axs:
                    ax.relim()
                    ax.autoscale_view()

                fig1.canvas.draw()
                
                with open(OUTPUT_FILE_PATH, 'w') as file:
                    file.write('Time (s), Thrust (N), CC Pressure (Pa), '+\
                        'mdot_out (kg/s), T_t (K), OF (dimless), k (dimless), mdot_in (kg/s)')

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

                ax10.clear()
                ax10.text(0.9,0.7, f'CC Pres: {engine_model.cc_pressure:.3f}', 
                            ha = 'right', transform=ax10.transAxes, fontsize = 'small')
                ax10.text(0.9,0.6, f'Choked?: {engine_model.choked}', 
                            ha = 'right', transform=ax10.transAxes, fontsize = 'small')

                if engine_model.pressure_build_time_evaluated:
                    ax10.text(0.9,0.9, f'Pres rise time: {engine_model.pressure_build_time:.3f}', 
                            ha = 'right', transform=ax10.transAxes, fontsize = 'small')
                    ax10.text(0.9,0.8, f'Peak thrust: {engine_model.thurst_at_pressure_peak:.3f}', 
                            ha = 'right', transform=ax10.transAxes, fontsize = 'small')
                    text_display_flag = True

                
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
                numerical_error_values.append(engine_model.numerical_error)

                step_error_values.append(step_error)
                step_time_values.append(step_size)

            


if __name__ == '__main__':
    if config["flags"]["sim_type"] == "solid":
        raise NotImplementedError("Solid simulation has not yet been implemented.")
    elif config["flags"]["sim_type"] == "hybrid":
        HybridBurnSimulator.sim_full_burn()
    elif config["flags"]["sim_type"] == "liquid":
        raise NotImplementedError("Liquid simulation has not yet been implemented.")
    else:
        raise RuntimeError("Invalid engine simulation type specified")
