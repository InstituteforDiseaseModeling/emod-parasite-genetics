import os
import numpy as np
from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from dtk.generic.climate import set_climate_constant

from simtools.SetupParser import SetupParser
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.ModBuilder import  ModFn, ModBuilder

from malaria.reports.MalariaReport import add_patient_report, add_survey_report, add_summary_report
from malaria.reports.MalariaReport import add_malaria_transmission_report

from dtk.interventions.input_EIR import add_InputEIR

from var_gene_interventions import add_var_gene_outbreak, add_malaria_outbreak

# General --------------------------------------------------------------------------------------------------------

years = 30 # length of simulation, in years
num_seeds = 1
num_infections = 10
num_hepatocytes = 1
report_start = 0
# exp_name = f'CoTransmission_no_seasonality_start_report_{report_start}'
exp_name = f'test_wide_immune_param_sweep'
# Setup ----------------------------------------------------------------------------------------------------------
# config_path = os.path.join('.', 'inputs','config.json')
# cb = DTKConfigBuilder.from_files(config_path)

cb = DTKConfigBuilder.from_defaults('MALARIA_SIM')

cb.update_params({
    "Demographics_Filenames": ["challenge_infection_demographics.json"],
    "Campaign_Filename": "campaign.json",
    "Malaria_Strain_Model": "FALCIPARUM_RANDOM_STRAIN",
    "Falciparum_MSP_Variants": 100,
    "Falciparum_Nonspecific_Types": 20,
    "Falciparum_PfEMP1_Variants": 1000,
    "Hyperimmunity_Halflife": 120,
    "Base_Sporozoite_Survival_Fraction": 0.1*num_hepatocytes,
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
    "Max_Individual_Infections": 100

})

set_climate_constant(cb)

#Experimental Design -------------------------------------------------------------------------------------------------
# def sweep_larval_habitat(cb, scale_factor) :
#     cb.update_params({"x_Temporary_Larval_Habitat": scale_factor })
#     return { 'larval_habitat_multiplier' : scale_factor}


# add_patient_report(cb)
def add_outbreaks(cb, n_infections, coverage, period):

    for i in np.arange(n_infections):

        add_malaria_outbreak(cb,start_days=[0], coverage = coverage, repetitions= years*12, tsteps_btwn_repetitions = period,
                             Ignore_Immunity=False)

    return {'Coverage': coverage, 'Period': period, 'Infection number': n_infections}


def set_antigen_space(cb, var_length, switch_rate,memory_level,nonspec_factor):

    cb.update_params({"Falciparum_PfEMP1_Variants": int(var_length),
                      "Antigen_Switch_Rate":switch_rate,
                      "Antibody_Memory_Level":memory_level,
                      "Nonspecific_Antigenicity_Factor":nonspec_factor})
    return {'PFalciparum_PfEMP1_Variants': int(var_length), 'Antigen_Switch_Rate': switch_rate, 'Antibody_Memory_Level': memory_level,
            'Nonspecific_Antigenicity_Factor':nonspec_factor}

infection_number_array = [1,2,5,10]
coverage_array = [1]
period_array = [180]
var_lengths = list(np.logspace(start=2, stop = 4, num= 10))
# switch_rates = list(np.logspace(start=-11, stop = -9, num= 10))
memory_levels = list(np.linspace(start=0.01,stop = 0.35, num = 10))
# nonspec_factors = list(np.linspace(start=0.01,stop = 0.5, num = 10))
hyper_immune_decay_rates = list(np.linspace(start=0.1,stop = 0.5, num = 10))

builder = ModBuilder.from_list(
    [
            [
                ModFn(DTKConfigBuilder.set_param,'Run_Number',seed),
                ModFn(add_outbreaks,n_infections, coverage, period),
                ModFn(set_antigen_space,var_length,switch_rate,memory_level,nonspec_factor)
            ]
        for seed in range(num_seeds)
        for n_infections in infection_number_array
        for coverage in coverage_array
        for period in period_array
        for var_length in var_lengths
        for switch_rate in switch_rates
        for memory_level in memory_levels
        for nonspec_factor in nonspec_factors
        for hyper_immune_decay_rate in hyper_immune_decay_rates

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
