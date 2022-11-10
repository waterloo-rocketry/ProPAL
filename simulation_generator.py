import csv
import os
from tank_blowdown_model import NOS_tank
from engine_model import CombustionChamberModel, N2O_HTPB_ThermochemistryModel

import math
import decimal
import numpy as np

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

    iters = []
    for _,val in params.items():
      temp = []
      # get precision 
      precision = self.get_precision(val[2])
      for x in np.linspace(val[0], val[1], int((val[1]-val[0])/val[2] + 1)):
        temp.append(round(x, precision))
      iters.append(temp)
      # writer.writerow(temp)

    for item in product(*iters):
      # result = cantera_function(*item)
      writer.writerow([*item])

    print("writing done")
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
  def __init__(self) -> None: 
    pass

sg = SimulationGenerator('.')

stuff = {
  'OF_ratio': [1.000, 5.000, .2],
  'temperature_K': [290, 310, 0.5], 
  'pressure_Pa': [10000, 20000, 100]
}

sg.generate_lookup(stuff, N2O_HTPB_ThermochemistryModel.sim_gas_mixture_combustion_temp)
# sg.remove_lookup(stuff)