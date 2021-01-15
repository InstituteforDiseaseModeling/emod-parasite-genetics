import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import seaborn as sns
import matplotlib as mpl
from simtools.Analysis.AnalyzeManager import AnalyzeManager
from simtools.Analysis.BaseAnalyzers import BaseAnalyzer
from simtools.SetupParser import SetupParser
from scipy import interpolate
import os
import json

#read json to df
def json_to_df(fp):
    with open(fp) as f:
        data = json.load(f)
    df = pd.DataFrame(data['next_point']['data'])
    return df

def plot_df_sns(df,param_names,output_path):
    for i,p in enumerate(param_names):
        splot = sns.scatterplot(data=df, x = df[p], y = df['Results'])
        fig = splot.get_figure()
        fig_name = f'recalibration_results_for_inf2_{p}'
        fig.savefig(os.path.join(output_path,fig_name))
        fig.clf()

# def pairplot_df():

if __name__ == '__main__' :

    output_path = os.path.join(os.path.expanduser('~'),'Dropbox (IDM)','Malaria Team Folder','projects',
                               'parasite_genetics','DTK','VarGenes','outputs','recalibration')
    params = [
        {
            'Name': 'Falciparum_Nonspecific_Types',
            'Dynamic': True,
            'Guess': 20,
            'Min': 1,
            'Max': 100
        },
        {
            'Name': 'Falciparum_PfEMP1_Variants',
            'Dynamic': True,
            'Guess': 1070,
            'Min': 1000,
            'Max': 5000
        },
        {
            'Name': 'Antibody_Memory_Level',
            'Dynamic': True,
            'Guess': 0.15,
            'Min': 0.1,
            'Max': 0.35
        },
        {
            'Name': 'Hyperimmunity_Halflife',
            'Dynamic': True,
            'Guess': 120,
            'Min': 60,
            'Max': 180
        }
    ]

    param_names = [p['Name'] for p in params]
    param_min = [p['Min'] for p in params]
    param_max = [p['Max'] for p in params]


    num_infs = [2,5,10]
    agg_df = pd.DataFrame()
    fig,axes = plt.subplots(1,4, figsize = (15,5))
    for num_inf in num_infs:

        fp = os.path.join(os.getcwd(),'..',f'ImmunityCalib_{num_inf}inf_20samps_25iters_test_dtype','iter21','IterationState.json')
        df = json_to_df(fp)
        df['num_inf'] = num_inf
        agg_df = pd.concat([agg_df,df], sort=True)
    for i, p in enumerate(param_names):
        splot = sns.scatterplot(ax = axes[i], data=agg_df, x=agg_df[p], y=agg_df['Results'], hue = agg_df['num_inf'],
                                alpha = 0.5, palette = ['#66BB6A','#00ACC1','#1A237E'])

    fig_name = 'recalibration_results_for_inf__num_sweep'
    fig.savefig(os.path.join(output_path, fig_name))
    plt.show()
    fig.clf()


    #make facet grid plot for all params
    # load in the relevant json
