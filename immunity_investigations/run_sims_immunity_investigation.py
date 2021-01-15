from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from dtk.generic.climate import set_climate_constant

from simtools.SetupParser import SetupParser
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.ModBuilder import  ModFn, ModBuilder

from immunity_investigations.var_gene_interventions import add_var_gene_outbreak
# General --------------------------------------------------------------------------------------------------------

years = 5 # length of simulation, in years
num_seeds = 5
report_start = 0
# exp_name = f'CoTransmission_no_seasonality_start_report_{report_start}'
exp_name = 'simultaneous_heterologous_challenge'
# Setup ----------------------------------------------------------------------------------------------------------
# config_path = os.path.join('.', 'inputs','config.json')
# cb = DTKConfigBuilder.from_files(config_path)

cb = DTKConfigBuilder.from_defaults('MALARIA_SIM')

cb.update_params({
    "Demographics_Filenames": ["challenge_infection_demographics.json"],
    "Campaign_Filename": "campaign.json",
    "Malaria_Strain_Model": "FALCIPARUM_FIXED_STRAIN",
    "Falciparum_MSP_Variants": 100,
    "Falciparum_Nonspecific_Types": 20,
    "Falciparum_PfEMP1_Variants": 1000,
    "PKPD_Model": "CONCENTRATION_VERSUS_TIME",
    "Number_Basestrains": 1,
    "Genome_Markers": [

    ],
    "Vector_Sampling_Type": "TRACK_ALL_VECTORS",
    "x_Base_Population":0.5,
    "Enable_Disease_Mortality":0,
    "Enable_Vital_Dynamics":0,
    "Simulation_Duration":(years*365)+1,
    "Enable_Malaria_CoTransmission": 0,
    "Max_Individual_Infections": 10

})

set_climate_constant(cb)

#Experimental Design -------------------------------------------------------------------------------------------------
# def sweep_larval_habitat(cb, scale_factor) :
#     cb.update_params({"x_Temporary_Larval_Habitat": scale_factor })
#     return { 'larval_habitat_multiplier' : scale_factor}

infection_a = {'irbc_type': [
                        1,74,147,220,293,366,439,512,585,658,731,804,877,950,23,96,169,242,315,388,461,534,
                        607,680,753,826,899,972,45,118,191,264,337,410,483,556,629,702,775,848,921,994,67,
                        140,213,286,359,432,505,578
                    ],
                'minor_epitope_type': [
                        2,1,1,4,3,0,2,1,0,2,3,3,1,3,4,1,4,
                        1,1,4,3,1,4,4,1,4,1,1,1,0,1,2,3,0,
                        2,0,0,2,2,2,3,2,3,1,4,4,1,0,1,2
                    ]
}
infection_b = {'irbc_type': [2, 75, 148, 221, 294, 367, 440, 513, 586, 659, 732, 805, 878, 951, 24, 97, 170, 243, 316, 389, 462, 535, 608, 681, 754, 827, 900, 973, 46, 119, 192, 265, 338, 411, 484, 557, 630, 703, 776, 849, 922, 995, 68, 141, 214, 287, 360, 433, 506, 579],
                'minor_epitope_type': [
                        2,1,1,4,3,0,2,1,0,2,3,3,1,3,4,1,4,
                        1,1,4,3,1,4,4,1,4,1,1,1,0,1,2,3,0,
                        2,0,0,2,2,2,3,2,3,1,4,4,1,0,1,2
                    ]
}

# add_patient_report(cb)
add_var_gene_outbreak(cb,start_days = [1],coverage=1, repetitions=1,tsteps_btwn_repetitions=30,nodeIDs=[1],
                      node_property_restrictions=None,
                      ind_property_restrictions=None,
                      irbc_type = infection_a['irbc_type'],
                      msp_type = 1,
                      minor_epitope_type = [
                        2,1,1,4,3,0,2,1,0,2,3,3,1,3,4,1,4,
                        1,1,4,3,1,4,4,1,4,1,1,1,0,1,2,3,0,
                        2,0,0,2,2,2,3,2,3,1,4,4,1,0,1,2
                    ]
)
add_var_gene_outbreak(cb,start_days = [1],coverage=1, repetitions=1,tsteps_btwn_repetitions=30,nodeIDs=[1],
                      node_property_restrictions=None,
                      ind_property_restrictions=None,
                      irbc_type = infection_b['irbc_type'],
                      msp_type = 1,
                      minor_epitope_type = [
                        2,1,1,4,3,0,2,1,0,2,3,3,1,3,4,1,4,
                        1,1,4,3,1,4,4,1,4,1,1,1,0,1,2,3,0,
                        2,0,0,2,2,2,3,2,3,1,4,4,1,0,1,2
                    ]
)


builder = ModBuilder.from_list(
    [
            [
                ModFn(DTKConfigBuilder.set_param,'Run_Number',seed)
            ]
        for seed in range(num_seeds)
    ]
)




# Run args
run_sim_args = {'config_builder': cb,
                'exp_name': exp_name,
                'exp_builder': builder
                }

if __name__ == "__main__":

    if not SetupParser.initialized:
        SetupParser.init('HPC')

    exp_manager = ExperimentManagerFactory.init()
    exp_manager.run_simulations(**run_sim_args)
    exp_manager.wait_for_finished(verbose=True)
    assert (exp_manager.succeeded())
