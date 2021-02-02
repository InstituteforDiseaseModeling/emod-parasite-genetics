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

mpl.rcParams['pdf.fonttype'] = 42
if not SetupParser.initialized:
    SetupParser.init('HPC')

wdir = os.path.join(os.getcwd(),'output')

class PfPR_Analyzer(BaseAnalyzer):
    def __init__(self, output_fname):
        super(PfPR_Analyzer, self).__init__()
        self.filenames = ['output/InsetChart.json']
        self.output_fname = output_fname
        self.channels = ['PCR Parasite Prevalence']


    def select_simulation_data(self, data, simulation):
        simdata = pd.DataFrame()

        year_to_report = 25

        for channel in self.channels:

            value = data[self.filenames[0]]['Channels'][channel]['Data'][year_to_report*365::]
            simdata[channel] = pd.Series(value)


        # simdata[self.tag] = simulation.tags[self.tag]

        simdata['id'] = simulation.id
        simdata['seed'] = simulation.tags['Run_Number']


        return simdata

    def finalize(self, all_data):

        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("No data have been returned... Exiting...")
            return

        cmap = cm.get_cmap('viridis',1).colors

        fig, axes = plt.subplots(ncols=1, nrows=len(self.channels))
        cmap = cm.get_cmap('viridis',len(selected)).colors

        for i,sim in enumerate(selected):
            for j,ch in enumerate(self.channels):
                axes.plot(sim[ch], color='b',alpha = 0.5)
                axes.set_ylabel(ch)
                axes.set_ylim(0,1)
                axes.set_xlabel('Day')

        fig_name = f'Garki_Immune_Variation_All_Random_last5'
        plt.tight_layout()
        # fig.suptitle('Seasonal scenario with IRS')
        plt.savefig(
            rf'C:\Users\jorussell\Dropbox (IDM)\Malaria Team Folder\projects\parasite_genetics\DTK\Epi_Validity_1to1_biting\figures\{fig_name}.eps')
        plt.savefig(
            rf'C:\Users\jorussell\Dropbox (IDM)\Malaria Team Folder\projects\parasite_genetics\DTK\Epi_Validity_1to1_biting\figures\{fig_name}.png')
        plt.show()


        # df = pd.concat(selected).reset_index(drop=True)
        # df.to_csv('%s_inset_chart.csv' % self.output_fname)

if __name__ == '__main__' :
    # import json
    # with open('C:\git\magude\pickup_realistic\inputs\Demographics\demo_exe_224.json') as json_file:
    #     datafile = json.load(json_file)
    #     print([x['NodeID'] for x in datafile['Nodes']])


    expids = ['c4c4ee27-1555-eb11-a2dd-c4346bcb7271']


    expnames = ['base_demo']
    channel_name = 'PCR Parasite Prevalence'
    for expname, expid in zip(expnames, expids) :
        output_fname = os.path.join(wdir, expid,expname)
        am = AnalyzeManager(exp_list = expids,
                            analyzers=PfPR_Analyzer(output_fname))
        am.analyze()
