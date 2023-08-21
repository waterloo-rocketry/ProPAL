import csv
import os
import sys
sys.path.insert(0,'.')
from ThermochemWrappers import N2O_HTPB_ThermochemistryModel

import math
import decimal
import numpy as np
from time import perf_counter

from itertools import product

#########################################################################################################

class SimulationGenerator:
  # a class that can generate any number of LUTs at the given path
  def __init__(self, folderpath = '.') -> None:
    # create folderpath if does not exist
    self.folderpath = folderpath
    if folderpath != '.' and not os.path.exists(folderpath):
      os.mkdir(folderpath)
    print('SimulationGenerator will save to:', os.path.abspath(self.folderpath))
    
  def get_precision(self, number):
    if number >= 1: return 0
    return -decimal.Decimal(str(number)).as_tuple().exponent

  # params: a dictionary of keys (str) and ranges ([min, max, stepsize] arrays)
  # cantera_function: a function that should take the same # of params as your params param does (lol). This is the
  #   fn you want to optimize away
  # filepath: the path to save the file to

  def generate_lookup(self, params, cantera_function, file_prefix = None):
    filename = self.folderpath + '/SG-' + file_prefix + '-'.join([key for key in params]) + '---' + '-'.join([str(v) for _,value in params.items() for v in value]) + '.csv'
    print(filename)
    if os.path.exists(filename):
      print("file exists at", filename)
      return 
    
    f = open(filename, 'w')
    writer = csv.writer(f)

    # Writing header of generated file
    writer.writerow([*[key for key in params], 'P', 'CV', 'CP', 'T', 'M_avg', 'rho'])

    # write all possible values for each key
    iters = []
    for _,val in params.items():
      temp = []
      # get precision 
      precision = self.get_precision(val[2])
      for x in np.arange(val[0], val[1]+val[2], val[2]):
        temp.append(round(x, precision))
      iters.append(temp)
    
    number_of_calculations = 1
    for sublist in iters:
      number_of_calculations *= len(sublist)
    
    print('Generating the target table will take ' + str(number_of_calculations) + ' calculations')
    print(iters)

    prev_value = iters[0][0]
    last_value = iters[0][-1]

    start_time = perf_counter()
    for item in product(*iters):
      if item[0] > prev_value:
        prev_value = item[0]
        print(prev_value, '/', last_value)
      result = cantera_function(*item)
      writer.writerow([*item, result.P, result.cv, result.cp, result.T, result.mean_molecular_weight, result.density])

    print("writing done")
    print('Cantera simulation generator time: ' + str(perf_counter() - start_time))
    f.close()

    
  # function below is solely for testing purposes
  def remove_lookup(self, params): 
    filename = self.folderpath + '/SG-' + '-'.join([key for key in params]) + '---' + '-'.join([str(v) for _,value in params.items() for v in value]) + '.csv'
    if os.path.exists(filename):
      os.remove(filename)
      print(filename, "removed.")
    else: 
      print(filename, "does not exist at", self.folderpath + '/')


class MockCanteraResult:
    def __init__(self, iP, icv, icp, iT, iM_avg, i_rho) -> None:
      self.P = iP
      self.cv = icv
      self.cp = icp
      self.T = iT
      self.mean_molecular_weight = iM_avg
      self.density = i_rho

class LookUpTable:

  def __init__(self, filename, step_sizes=None) -> None: 
    f = open(filename, 'r')
    self.reader = csv.reader(f)


    self.columns = next(self.reader)
    self.num_columns = len(self.columns)

    self.num_input_params = filename.split('-').index('') - 1
    print(self.num_input_params)
    



    self.table = dict()
    for row in self.reader:
      # print(len(row))
      if len(row) == 0:
        continue


      self.table[tuple([float(row[i]) for i in range(self.num_input_params)])] = np.array([float(row[i]) for i in range(self.num_input_params, len(row))])

    # self.table = np.array(temp_table)
    # print(self.table[(1.0,290.0,10000)])
    f.close()

    if step_sizes:
      self.step_sizes = step_sizes
    else:
      self.infer_step_sizes_and_starts()

  def infer_step_sizes_and_starts(self):

    #extract all key values
    all_p1_keys = [key[0] for key in self.table.keys()]
    all_p2_keys = [key[1] for key in self.table.keys()]
    all_p3_keys = [key[2] for key in self.table.keys()]

    # print(all_p1_keys)

    # Extract only unique values and sort them
    unique_p1_keys_sorted = list(set(all_p1_keys))
    unique_p2_keys_sorted = list(set(all_p2_keys))
    unique_p3_keys_sorted = list(set(all_p3_keys))


    unique_p1_keys_sorted.sort()
    unique_p2_keys_sorted.sort()
    unique_p3_keys_sorted.sort()

    # print(unique_p1_keys_sorted)

    # infer step size from sorted arrays
    self.param1_step = round(unique_p1_keys_sorted[1] - unique_p1_keys_sorted[0], 5)
    self.param2_step = round(unique_p2_keys_sorted[1] - unique_p2_keys_sorted[0], 5)
    self.param3_step = round(unique_p3_keys_sorted[1] - unique_p3_keys_sorted[0], 5)

    self.param1_start = unique_p1_keys_sorted[0]
    self.param2_start = unique_p2_keys_sorted[0]
    self.param3_start = unique_p3_keys_sorted[0]


  # def lookup_old(self, params_to_search):
  #   # perform trilinear interpolation
  #   if self.table[params_to_search]:
  #     return self.table[params_to_search]
  #   # else 
  #   nearest_points = [ [] for _ in params_to_search ]
  #   # how do we find nearest points without searching the entire dictionary??
  #   # aaaaaaaaaaa
  #   # Solution: could be something calculated from step size and starting point

  #   iters = [] # declared the same as in the above fn
  #   # for _,val in params.items():
  #   #   temp = []
  #   #   # get precision 
  #   #   precision = self.get_precision(val[2])
  #   #   for x in np.arange(val[0], val[1]+val[2], val[2]):
  #   #     temp.append(round(x, precision))
  #   #   iters.append(temp)

  #   index = 0
  #   for key, val in params_to_search:
  #     for j in len(iters[index]):
  #       if iters[index][j] < params_to_search[val] and iters[index][j] + step_size > params_to_search[val]:
  #         nearest_points[index].append(iters[index][j])
  #       elif iters[index][j] > params_to_search[val] and iters[index][j] - step_size < params_to_search[val]:
  #         nearest_points[index].append(iters[index][j])
  #     index += 1
      
  #   pass

  def lookup(self, p1,p2,p3) -> tuple:

    # Manual implementation of trilinear interpolation algorithm. 

    p1_val_lower = round(math.floor((p1 - self.param1_start) /self.param1_step) * self.param1_step, 5) + self.param1_start
    p1_val_upper = round(math.ceil((p1 - self.param1_start)/self.param1_step) * self.param1_step, 5) + self.param1_start

    p2_val_lower = round(math.floor((p2 - self.param2_start)/self.param2_step) * self.param2_step, 5) + self.param2_start
    p2_val_upper = round(math.ceil((p2 - self.param2_start)/self.param2_step) * self.param2_step + self.param2_start, 5)

    p3_val_lower = round(math.floor((p3 - self.param3_start)/self.param3_step) * self.param3_step, 5) + self.param3_start
    p3_val_upper = round(math.ceil((p3 - self.param3_start)/self.param3_step) * self.param3_step, 5) + self.param3_start

    p1_difference = (p1 - p1_val_lower)/(p1_val_upper - p1_val_lower)
    p2_difference = (p2 - p2_val_lower)/(p2_val_upper - p2_val_lower)
    p3_difference = (p3 - p3_val_lower)/(p3_val_upper - p3_val_lower)

    C_000 = self.table[p1_val_lower, p2_val_lower, p3_val_lower]

    C_100 = self.table[p1_val_upper, p2_val_lower, p3_val_lower]
    C_010 = self.table[p1_val_lower, p2_val_upper, p3_val_lower]
    C_001 = self.table[p1_val_lower, p2_val_lower, p3_val_upper]

    C_110 = self.table[p1_val_upper, p2_val_upper, p3_val_lower]
    C_011 = self.table[p1_val_lower, p2_val_upper, p3_val_upper]
    C_101 = self.table[p1_val_upper, p2_val_lower, p3_val_upper]

    C_111 = self.table[p1_val_upper, p2_val_upper, p3_val_upper]

    C_00 = self.linear_interpolate_structs(C_000, C_100, p1_difference)
    C_10 = self.linear_interpolate_structs(C_001, C_101, p1_difference)
    C_01 = self.linear_interpolate_structs(C_010, C_110, p1_difference)
    C_11 = self.linear_interpolate_structs(C_011, C_111, p1_difference)

    C_0 = self.linear_interpolate_structs(C_00, C_01, p2_difference)
    C_1 = self.linear_interpolate_structs(C_10, C_11, p2_difference)

    C_final = self.linear_interpolate_structs(C_0, C_1, p3_difference)

    return MockCanteraResult(*C_final)


  def linear_interpolate_structs(self, tuple_0, tuple_1, argument):

    if len(tuple_1) != len(tuple_0):
      raise RuntimeError('Dimension Mismatch in custom linear interpolation function')
    

    ret_list = [0 for _ in tuple_0]

    for idx in range(len(ret_list)):
      ret_list[idx] = tuple_0[idx] * (1 - argument) + tuple_1[idx]*argument
    
    return ret_list


class Sim:
  def __init__(self, lookup_table_path=None) -> None:
    if lookup_table_path == None:
      self.generate_from_lookup = False
    else:
      self.generate_from_lookup = True
    pass
  


def dev_case_1():
  sg = SimulationGenerator('.')

  stuff = {
    'OF_ratio': [1.000, 5.000, .3],
    'temperature_K': [290, 310, 0.5], 
    'pressure_Pa': [10000, 20000, 250]
  }

  # sg.generate_lookup(stuff, N2O_HTPB_ThermocheNmistryModel.sim_gas_mixture_combustion_temp)
  # sg.remove_lookup(stuff)

  lut = LookUpTable("OptimizationTemp/SG-OF_ratio-temperature_K-pressure_Pa---1.0-5.0-0.3-290-310-0.5-10000-20000-100.csv", [0.3, 0.5, 100])


def dev_case_2():
  sg = SimulationGenerator('.')

  table_params = {
    'OF_ratio': [1.000, 5.000, .3],
    'temperature_K': [290, 310, 0.5], 
    'pressure_Pa': [10000, 20000, 250]
  }

  sg.generate_lookup(table_params, N2O_HTPB_ThermochemistryModel.sim_gas_mixture_combustion_temp, file_prefix='TargetedTable')

def dev_case_3():
  sg = SimulationGenerator('.')

  table_params = {
    'OF_ratio': [0.1, 15, 0.25],
    'temperature_K': [273.3, 273.3, 0.1], 
    'pressure_Pa': [1, 1e7, 25000]
  }

  sg.generate_lookup(table_params, N2O_HTPB_ThermochemistryModel.sim_gas_mixture_combustion_temp, file_prefix='TargetedTable')
  
def dev_case_2_teardown():

  sg = SimulationGenerator('.')

  table_params = {
    'OF_ratio': [1.000, 5.000, .3],
    'temperature_K': [290, 310, 0.5], 
    'pressure_Pa': [10000, 20000, 250]
  }

  sg.remove_lookup(table_params)

def dev_case_4():


  lut_creation_start = perf_counter()
  lut = LookUpTable("OptimizationTemp/SG-TargetedTableOF_ratio-temperature_K-pressure_Pa---0.1-15-0.25-273.3-273.3-0.1-1-10000000.0-25000.csv")

  print("Creation time for lookup table: " + str(perf_counter() - lut_creation_start))

  result = lut.lookup(3, 273.35, 101325)
  # print(lut.table)

  print("T: " + str(result.T))
  print("P: " + str(result.P))
  print("cv: " + str(result.cv))
  print("cp: " + str(result.cp))
  print("M_mean: " + str(result.mean_molecular_mass))
  print("rho: " + str(result.density))


def dev_case_5():
  sg = SimulationGenerator('.')

  table_params = {
    'OF_ratio': [0.05, 15, 0.05],
    'temperature_K': [273.3, 273.3, 0.1], 
    'pressure_Pa': [1, 1e7, 2500]
  }

  sg.generate_lookup(table_params, N2O_HTPB_ThermochemistryModel.sim_gas_mixture_combustion_temp, file_prefix='IncreasedResolution')


if __name__ == "__main__":
  dev_case_4()
