import csv
import os
from tank_blowdown_model import NOS_tank
from engine_model import CombustionChamberModel

import math
import numpy as np

#########################################################################################################

class SimulationGenerator:
  # a class that can generate any number of LUTs at the given path
  def __init__(self, folderpath = '.') -> None:
    # create folderpath if does not exist
    self.folderpath = folderpath
    if folderpath != '.' and not os.path.exists(folderpath):
      os.mkdir(folderpath)
    print('SimulationGenerator will save to:', os.path.abspath(self.folderpath))
    
  
  # params: a dictionary of keys (str) and ranges ([min, max, stepsize] arrays)
  # filepath: the path to save the file to

  def generate_lookup(self, params):
    self.filename = self.folderpath + '/SG-' + '-'.join([key for key in params]) + '---' + '-'.join([str(v) for _,value in params.items() for v in value]) + '.csv'
    print(self.filename)
    if os.path.exists(self.filename):
      print("file exists at", self.filename)
      return 
    
    f = open(self.filename, 'w')
    writer = csv.writer(f)
    writer.writerow([key for key in params])

    for key,val in params.items():
      for x in np.arange(val[0], val[1]+val[2], val[2]):
        print(x)
    
    
    


class LookUpTable:
  def __init__(self) -> None: 
    pass

# k imma make it a CSV actually bc import to excel and also lots of values r easier to parse with it
# NOS_tank_params = [0.04, 55, 288, 40, 0.15]
# SIM_STEP = 0.025
# tank_model = NOS_tank(*NOS_tank_params)
# comb_chamber_model = CombustionChamberModel(D_0=0.05, L=0.6, _fuel_density=1100)

# ox_massflow = tank_model.massflow
# comb_chamber_model.sim_comubstion(oxidizer_massflow=ox_massflow, delta_time=SIM_STEP)
# fuel_massflow = comb_chamber_model.get_fuel_massflow()
# # OF_ratio = ox_massflow/fuel_massflow

# bar_to_Pa = 100000

# pressure_Pa = bar_to_Pa*(tank_model.pressure/2)
sg = SimulationGenerator('.')

stuff = {
  'OF_ratio': [1.000, 5.000, 0.002],
  'temperature_K': [290, 310, 0.5], 
  # 'pressure_Pa': [10000, 20000, 0.2]
}
sg.generate_lookup(stuff)
