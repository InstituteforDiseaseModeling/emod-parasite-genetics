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


def smooth(x, window_len=10, window='hanning'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len < 3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    s = np.r_[x[window_len - 1:0:-1], x, x[-2:-window_len - 1:-1]]
    # print(len(s))
    if window == 'flat':  # moving average
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.' + window + '(window_len)')

    y = np.convolve(w / w.sum(), s, mode='valid')
    return y

class PfPR_Analyzer(BaseAnalyzer):
    def __init__(self, output_fname):
        super(PfPR_Analyzer, self).__init__()
        self.filenames = ['output/InsetChart.json']
        self.output_fname = output_fname
        self.channels = ['PCR Parasite Prevalence']


    def select_simulation_data(self, data, simulation):
        simdata = pd.DataFrame()

        year_to_report = 27

        for channel in self.channels:

            value = data[self.filenames[0]]['Channels'][channel]['Data'][year_to_report*365::]
            simdata[channel] = pd.Series(value)


        # simdata[self.tag] = simulation.tags[self.tag]
        # simdata['annual EIR'] = [x*365 for x in simdata['Daily EIR']]
        simdata['id'] = simulation.id
        # simdata['LHM'] = simulation.tags['larval_habitat_multiplier']
        simdata['seed'] = simulation.tags['Run_Number']
        simdata['period'] = simulation.tags['Period']
        simdata['infection_number'] = simulation.tags['Infection number']

        return simdata

    def finalize(self, all_data):

        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("No data have been returned... Exiting...")
            return

        cmap = cm.get_cmap('viridis',10).colors

        fig, axes = plt.subplots(ncols=1, nrows=len(self.channels))
        cmap = cm.get_cmap('viridis',len(selected)).colors

        for i,sim in enumerate(selected):
            for j,ch in enumerate(self.channels):
                axes.plot(smooth(sim[ch]), color=cmap[i],alpha = 0.5)
                axes.set_ylabel(ch)
        fig_name = f'base_results'
        # plt.xlabel('day')
        # axes.set_ylabel('PfPR')
        # ax.set_ylim(0, 1)

        plt.tight_layout()
        # fig.suptitle('Seasonal scenario with IRS')
        plt.savefig(
            rf'C:\Users\jorussell\Dropbox (IDM)\Malaria Team Folder\projects\parasite_genetics\DTK\VarGenes\outputs\memory_sweep\{fig_name}.eps')
        plt.savefig(
            rf'C:\Users\jorussell\Dropbox (IDM)\Malaria Team Folder\projects\parasite_genetics\DTK\VarGenes\outputs\memory_sweep\{fig_name}.png')
        plt.show()


        # df = pd.concat(selected).reset_index(drop=True)
        # df.to_csv('%s_inset_chart.csv' % self.output_fname)

if __name__ == '__main__' :
    # import json
    # with open('C:\git\magude\pickup_realistic\inputs\Demographics\demo_exe_224.json') as json_file:
    #     datafile = json.load(json_file)
    #     print([x['NodeID'] for x in datafile['Nodes']])


    expids = ['eea8e1c8-8218-eb11-a2c7-c4346bcb1553']

    simids = ['f1a8e1c8-8218-eb11-a2c7-c4346bcb1553']

    expnames = ['base_demo']
    channel_name = 'PCR Parasite Prevalence'
    for expname, simid in zip(expnames, simids) :
        output_fname = os.path.join(wdir, simid,expname)
        am = AnalyzeManager(sim_list = simids,
                            analyzers=PfPR_Analyzer(output_fname))
        am.analyze()
