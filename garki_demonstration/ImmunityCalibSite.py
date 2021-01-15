import logging
from immunity_investigations.Helpers import reference_prevalence

from immunity_investigations.analyzers.memory_sweep_analyzer import MemorySweepAnalyzer
from abc import ABCMeta
from calibtool.CalibSite import CalibSite

from malaria.study_sites.site_setup_functions import config_setup_fn, update_params, add_var_outbreak_fn


logger = logging.getLogger(__name__)

class ImmunityCalibSite(CalibSite):
    """
    An abstract class that implements a simulation setup fro comparison of baseline epi/prevalence behavior under different assumptions of immune model parameterization:
    - Sweeps of:
        - "Falciparum_PfEMP1_variants"
        - "Antibody_Memory_Level"
        - "Hyperimmunity_Halflife"
    """

    __metaclass__ = ABCMeta

    def __init__(self,reference_fname, run_duration, analyze_duration,start_analysis,n_infections,**kwargs):


        self.metadata = {
            'fname': reference_fname,
            'analyze_duration': 180,
            'run_duration': 365*run_duration,
            'start_date': 9000,
            'n_infections': n_infections
        }

        super(ImmunityCalibSite, self).__init__(self.metadata['n_infections'])


    def get_reference_data(self, reference_type):

        reference_data = reference_prevalence(self.metadata)

        return reference_data

    def get_setup_functions(self):
        setup_fns = [config_setup_fn(duration=self.metadata['run_duration']),
                     update_params({'Demographics_Filenames': ['challenge_infection_demographics.json'],
                                    'Climate_Model': 'CLIMATE_CONSTANT',
                                    "Campaign_Filename": "campaign.json",
                                    "Malaria_Strain_Model": "FALCIPARUM_RANDOM_STRAIN",
                                    "Falciparum_MSP_Variants": 100,
                                    "Falciparum_Nonspecific_Types": 20,
                                    "Falciparum_PfEMP1_Variants": 1000,
                                    "Hyperimmunity_Halflife": 120,
                                    "PKPD_Model": "CONCENTRATION_VERSUS_TIME",
                                    "Number_Basestrains": 1,
                                    "Genome_Markers": [
                                    ],
                                    "Vector_Sampling_Type": "TRACK_ALL_VECTORS",
                                    "x_Base_Population":0.5,
                                    "Enable_Disease_Mortality":0,
                                    "Enable_Vital_Dynamics":0,
                                    "Enable_Malaria_CoTransmission": 0,
                                    "Max_Individual_Infections": 100
                                    }),
                     add_var_outbreak_fn(start_day = 0,n_infections = self.metadata['n_infections'], coverage=1, period = 180,
                                         repetitions=(self.metadata['run_duration']*2/365))

                     ]

        return setup_fns

    def get_analyzers(self):
        return [
            MemorySweepAnalyzer(site=self, analyze_duration=self.metadata['analyze_duration'], start_date = self.metadata['start_date'])]
