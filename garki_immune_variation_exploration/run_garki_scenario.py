from dtk.utils.core.DTKConfigBuilder import DTKConfigBuilder
from dtk.generic.climate import set_climate_constant
from dtk.vector.species import update_species_param, set_species_param
from simtools.SetupParser import SetupParser
from simtools.ExperimentManager.ExperimentManagerFactory import ExperimentManagerFactory
from simtools.ModBuilder import  ModFn, ModBuilder
from malaria.reports.MalariaReport import add_summary_report
from immunity_investigations.var_gene_interventions import add_var_gene_outbreak
# General --------------------------------------------------------------------------------------------------------

years = 30 # length of simulation, in years
num_seeds = 100
report_year = 25
years_to_report = 1
# exp_name = f'CoTransmission_no_seasonality_start_report_{report_start}'
exp_name = 'FPG_epi_validity_test_ALL_RANDOM_100seeds'
# Setup ----------------------------------------------------------------------------------------------------------
# config_path = os.path.join('.', 'inputs','config.json')
cb = DTKConfigBuilder.from_files(config_name ='C:\git\emod-parasite-genetics\garki_immune_variation_exploration\inputs\config.json', campaign_name="C:\git\emod-parasite-genetics\garki_immune_variation_exploration\inputs\campaign_ParasiteGenetics.json")

# cb = DTKConfigBuilder.from_defaults('MALARIA_SIM')

cb.update_params({
    "Demographics_Filenames": ["Garki_single_demographics.json"],
    "Campaign_Filename": "campaign_ParasiteGenetics.json",
    "Malaria_Strain_Model": "FALCIPARUM_RANDOM_STRAIN",
    "Falciparum_MSP_Variants": 100,
    "Falciparum_Nonspecific_Types": 20,
    "Falciparum_PfEMP1_Variants": 1000,
    "PKPD_Model": "CONCENTRATION_VERSUS_TIME",
    "Number_Basestrains": 1,
    "Genome_Markers": [

    ],
    "Parasite_Genetics": {
        "Enable_FPG_Similarity_To_Base": 0,
        "Sporozoite_Life_Expectancy": 25,
        "Num_Sporozoites_In_Bite_Fail": 12,
        "Probability_Sporozoite_In_Bite_Fails": 0.5,
        "Num_Oocyst_From_Bite_Fail": 3,
        "Probability_Oocyst_From_Bite_Fails": 0.5,
        "Sporozoites_Per_Oocyst_Distribution": "GAUSSIAN_DISTRIBUTION",
        "Sporozoites_Per_Oocyst_Gaussian_Mean": 10000,
        "Sporozoites_Per_Oocyst_Gaussian_Std_Dev": 1000,
        "Crossover_Gamma_K": 2000.0,
        "Crossover_Gamma_Theta": 100.0,
        "Barcode_Genome_Locations": [
            311500,
            1116500,
            2140000,
            3290000,
            4323333, 4756667,
            5656667, 6123333,
            7056667, 7523333,
            8423333, 8856667,
            9790000, 10290000,
            11356667, 11923333,
            13156667, 13823333,
            15256667, 16023333,
            17690000, 18590000,
            20590000, 21690000
        ],
        "MSP_Genome_Location": 200000,
        "PfEMP1_Variants_Genome_Locations": [
            214333, 428667,
            958667, 1274333,
            1864900, 2139900, 2414900,
            2989900, 3289900, 3589900,
            4150000, 4410000, 4670000, 4930000,
            5470000, 5750000, 6030000, 6310000,
            6870000, 7150000, 7430000, 7710000,
            8250000, 8510000, 8770000, 9030000,
            9590000, 9890000, 10190000, 10490000,
            11130000, 11470000, 11810000, 12150000,
            12890000, 13290000, 13690000, 14090000,
            14950000, 15410000, 15870000, 16330000,
            17330000, 17870000, 18410000, 18950000,
            20150000, 20810000, 21470000, 22130000
        ],
        "Var_Gene_Randomness_Type": "ALL_RANDOM",
        "Neighborhood_Size_MSP": 4,
        "Neighborhood_Size_PfEMP1": 10
         },
    "Egg_Saturation_At_Oviposition": "SATURATION_AT_OVIPOSITION",
    "Malaria_Model": "MALARIA_MECHANISTIC_MODEL_WITH_PARASITE_GENETICS",
    "Inset_Chart_Reporting_Include_30Day_Avg_Infection_Duration": 0,
    "Serialized_Population_Reading_Type": "NONE",
    "Serialized_Population_Writing_Type": "NONE",
    "Vector_Sampling_Type": "TRACK_ALL_VECTORS",
    "Vector_Species_Params": [
      {
        "Acquire_Modifier": 0.8,
        "Adult_Life_Expectancy": 20,
        "Anthropophily": 0.65,
        "Aquatic_Arrhenius_1": 84200000000,
        "Aquatic_Arrhenius_2": 8328,
        "Aquatic_Mortality_Rate": 0.1,
        "Cycle_Arrhenius_1": 0,
        "Cycle_Arrhenius_2": 0,
        "Cycle_Arrhenius_Reduction_Factor": 0,
        "Days_Between_Feeds": 3,
        "Egg_Batch_Size": 100,
        "Genes": [],
        "Immature_Duration": 2,
        "Indoor_Feeding_Fraction": 0.5,
        "Infected_Arrhenius_1": 117000000000,
        "Infected_Arrhenius_2": 8336,
        "Infected_Egg_Batch_Factor": 0.8,
        "Infectious_Human_Feed_Mortality_Factor": 1.5,
        "Larval_Habitat_Types": {
          "LINEAR_SPLINE": {
            "Capacity_Distribution_Number_Of_Years": 1,
            "Capacity_Distribution_Over_Time": {
              "Times": [
                0,
                30.417,
                60.833,
                91.25,
                121.667,
                152.083,
                182.5,
                212.917,
                243.333,
                273.75,
                304.167,
                334.583
              ],
              "Values": [

              ]
            },
            "Max_Larval_Capacity": 316227766.01683795
          }
        },
        "Male_Life_Expectancy": 10,
        "Name": "gambiae",
        "Transmission_Rate": 0.9,
        "Vector_Sugar_Feeding_Frequency": "VECTOR_SUGAR_FEEDING_EVERY_DAY"
      }],
    "x_Base_Population": 0.5,
    "Base_Population_Scale_Factor":0.5,
    "Enable_Disease_Mortality":0,
    "Enable_Initial_Prevalence": 0,
    "Enable_FPG_Similarity_To_Base": 0,
    "Enable_Vital_Dynamics":1,
    "Sporozoite_Decay_Rate": 20,
    "Simulation_Duration":(years*365)+1,
    "Enable_Malaria_CoTransmission": 0,
    "Max_Individual_Infections": 10,

})

set_climate_constant(cb)

def update_vector_params(cb):
    cb.update_params({"Vector_Species_Names": ['gambiae']})
    set_species_param(cb, 'gambiae', 'Larval_Habitat_Types',
                      {"LINEAR_SPLINE": {
                          "Capacity_Distribution_Over_Time": {
                              "Times": [0.0, 30.417, 60.833, 91.25, 121.667, 152.083,
                                        182.5, 212.917, 243.333, 273.75, 304.167, 334.583],
                              "Values": [3, 0.8,  1.25, 0.1,  2.7, 8, 4, 35, 6.8,  6.5, 2.6, 2.1]
                          },
                          "Capacity_Distribution_Number_Of_Years": 1,
                          "Max_Larval_Capacity": pow(10, 8)
                      }})
    set_species_param(cb, "gambiae", "Indoor_Feeding_Fraction", 0.9)
    set_species_param(cb, "gambiae", "Adult_Life_Expectancy", 20)


#Experimental Design -------------------------------------------------------------------------------------------------
# def sweep_larval_habitat(cb, scale_factor) :
#     cb.update_params({"x_Temporary_Larval_Habitat": scale_factor })
#     return { 'larval_habitat_multiplier' : scale_factor}

#Reporting -----------------------------------------------------------------------------------------------------------

add_summary_report(cb,
                   start = report_year*365,
                   description='Annual_Report',
                   interval = 365,
                   nreports = years_to_report,
                   age_bins = [1, 4, 8, 18, 28, 43, 125],
                   parasitemia_bins = [0, 50, 200, 500, 2000000]
                   )
# add_var_gene_outbreak(cb,start_days = [1],coverage=1, repetitions=1,tsteps_btwn_repetitions=30,nodeIDs=[1],
#                       node_property_restrictions=None,
#                       ind_property_restrictions=None,
#                       irbc_type = [
#                         1,74,147,220,293,366,439,512,585,658,731,804,877,950,23,96,169,242,315,388,461,534,
#                         607,680,753,826,899,972,45,118,191,264,337,410,483,556,629,702,775,848,921,994,67,
#                         140,213,286,359,432,505,578
#                     ],
#                       msp_type = 1,
#                       minor_epitope_type = [
#                         2,1,1,4,3,0,2,1,0,2,3,3,1,3,4,1,4,
#                         1,1,4,3,1,4,4,1,4,1,1,1,0,1,2,3,0,
#                         2,0,0,2,2,2,3,2,3,1,4,4,1,0,1,2
#                     ]
# )

update_vector_params(cb)
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
