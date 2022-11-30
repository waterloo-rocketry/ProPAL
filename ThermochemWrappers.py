import cantera as ct

class N2O_HTPB_ThermochemistryModel:

    def __init__(self) -> None:
        self.gas_model = ct.Solution('Mevel2015-rocketry_modified.yaml')
    
    
    def sim_gas_mixture_combustion_temp(self, OF_ratio, temperature_K, pressure_Pa) -> float:

        # This requires a download of a NASA file and some alteration to
        # get it to compile correctly


        self.gas_model.TPY = temperature_K, pressure_Pa, f'N2O:{OF_ratio}, C4H6:1'

        self.gas_model.equilibrate('HP')

        return self.gas_model