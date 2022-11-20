import csv
import os
from tank_blowdown_model import NOS_tank
from engine_model import CombustionChamberModel, N2O_HTPB_ThermochemistryModel

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

  def generate_lookup(self, params, cantera_function):
    filename = self.folderpath + '/SG-' + '-'.join([key for key in params]) + '---' + '-'.join([str(v) for _,value in params.items() for v in value]) + '.csv'
    print(filename)
    if os.path.exists(filename):
      print("file exists at", filename)
      return 
    
    f = open(filename, 'w')
    writer = csv.writer(f)
    writer.writerow([*[key for key in params], 'P', 'CV', 'CP', 'T'])

    # write all possible values for each key
    iters = []
    for _,val in params.items():
      temp = []
      # get precision 
      precision = self.get_precision(val[2])
      for x in np.arange(val[0], val[1]+val[2], val[2]):
        temp.append(round(x, precision))
      iters.append(temp)

    prev_value = iters[0][0]
    last_value = iters[0][-1]

    start_time = perf_counter()
    for item in product(*iters):
      if item[0] > prev_value:
        prev_value = item[0]
        print(prev_value, '/', last_value)
      result = cantera_function(*item)
      writer.writerow([*item, result.P, result.cv, result.cp, result.T])

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
    


class LookUpTable:
  def __init__(self, filename, step_sizes) -> None: 
    f = open(filename, 'r')
    self.reader = csv.reader(f)
    self.step_sizes = step_sizes
    self.columns = next(self.reader)
    self.num_columns = len(self.columns)

    self.num_input_params = filename.split('-').index('') - 1
    print(self.num_input_params)

    self.table = dict()
    for row in self.reader:
      self.table[tuple([float(row[i]) for i in range(self.num_input_params)])] = np.array([float(row[i]) for i in range(self.num_input_params, len(row))])

    # self.table = np.array(temp_table)
    print(self.table[(1.0,290.0,10000)])
    f.close()

  def lookup(self, params_to_search):
    # perform trilinear interpolation
    if self.table[params_to_search]:
      return self.table[params_to_search]
    # else 
    nearest_points = [ [] for _ in params_to_search ]
    # how do we find nearest points without searching the entire dictionary??
    # aaaaaaaaaaa
    # Solution: could be something calculated from step size and starting point

    iters = [] # declared the same as in the above fn
    # for _,val in params.items():
    #   temp = []
    #   # get precision 
    #   precision = self.get_precision(val[2])
    #   for x in np.arange(val[0], val[1]+val[2], val[2]):
    #     temp.append(round(x, precision))
    #   iters.append(temp)

    index = 0
    for key, val in params_to_search:
      for j in len(iters[index]):
        if iters[index][j] < params_to_search[val] and iters[index][j] + step_size > params_to_search[val]:
          nearest_points[index].append(iters[index][j])
        elif iters[index][j] > params_to_search[val] and iters[index][j] - step_size < params_to_search[val]:
          nearest_points[index].append(iters[index][j])
      index += 1
      
    pass

class Sim:
  def __init__(self, lookup_table_path=None) -> None:
    if lookup_table_path == None:
      self.generate_from_lookup = False
    else:
      self.generate_from_lookup = True
    pass
  


sg = SimulationGenerator('.')

stuff = {
  'OF_ratio': [1.000, 5.000, .3],
  'temperature_K': [290, 310, 0.5], 
  'pressure_Pa': [10000, 20000, 250]
}

# sg.generate_lookup(stuff, N2O_HTPB_ThermochemistryModel.sim_gas_mixture_combustion_temp)
# sg.remove_lookup(stuff)

lut = LookUpTable("./SG-OF_ratio-temperature_K-pressure_Pa---1.0-5.0-0.3-290-310-0.5-10000-20000-100.csv", [0.3, 0.5, 100])