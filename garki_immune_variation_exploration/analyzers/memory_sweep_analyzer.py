
import os
import logging
import pandas as pd
import numpy as np
from scipy import interpolate
from numpy.linalg import norm

import matplotlib.pyplot as plt
from matplotlib import cm
import seaborn as sns

from calibtool import LL_calculators
from simtools.Analysis.BaseAnalyzers import BaseCalibrationAnalyzer
from simtools.Analysis.AnalyzeManager import AnalyzeManager
from simtools.Analysis.BaseAnalyzers import BaseAnalyzer
from simtools.SetupParser import SetupParser


logger = logging.getLogger(__name__)

wdir = os.path.join(os.getcwd(),'output')

class MemorySweepAnalyzer(BaseCalibrationAnalyzer):
    """
    Base class implementation for comparison of population prevalence profiles from modified immune model to base behavior
    """

    def __init__(self, site, weight=1, compare_fn=LL_calculators.euclidean_distance_pandas, analyze_duration = 1*365,start_date = 25*365):
        super().__init__(reference_data=site.get_reference_data('true prevalence'),
                         weight=weight,
                         filenames=['output/InsetChart.json'])

        self.channel = 'True Prevalence'
        self.compare_fn = compare_fn
        self.site_name = site.name
        self.duration = analyze_duration
        self.start_date = start_date

    def select_simulation_data(self, data, simulation):
        """
        Extract data from output data and accumulate in same bins as reference.
        """

        # Load data from simulation
        simdata = {}
        for file in self.filenames:
            simdata = {self.channel:

                               data[file]['Channels']['True Prevalence']['Data'][self.start_date:(self.start_date+self.duration)]
                       }

            simdata = pd.DataFrame(simdata)
            simdata = simdata.rename(columns={self.channel: 'True Prevalence'
                                              })
        return simdata

    @staticmethod
    def join_reference(sim, ref):
        sim.columns = sim.columns.droplevel(level = [0,1])  # drop sim 'sample' to match ref levels
        return pd.concat({'sim': sim, 'ref': ref}, axis=1).dropna()

    def compare(self, sample):
        """
        Assess the result per sample, in this case the likelihood
        comparison between simulation and reference data.
        """
        return self.compare_fn(self.join_reference(sample, self.reference_data))

    def finalize(self, all_data):
        """
        Calculate the output result for each sample.
        """
        selected = list(all_data.values())

        # Stack selected_data from each parser, adding unique (sim_id) and shared (sample) levels to MultiIndex
        combine_levels = ['sample', 'sim_id','channel']
        combined = pd.concat(selected, axis=1,
                             keys=[(s.tags.get('__sample_index__'), s.id) for s in all_data.keys()],
                             names=combine_levels)

        compare_results = combined.groupby(level='sample', axis=1).apply(self.compare)

        ref = self.reference_data.reset_index()


        head, tail = os.path.split(self.working_dir)
        iteration = int(tail.split('r')[-1])
        sample_index = [s.tags.get('__sample_index__') for s in all_data.keys()]

        for n in range(len(selected)):
            fig = plt.figure()
            ax = fig.gca()
            plot_df = selected[n].reset_index()

            ax.plot(plot_df['True Prevalence'], '--', color='b',
                    label='iter %d sample %d' % (iteration, sample_index[n]))
            ax.plot(ref['True_PfPR'],
                    '--', color='r', label='reference')
            ax.set_xlabel('Time')
            ax.set_ylabel('True PfPR')
            # ax.set_ylim(-0.02, max(max(plot_df['fraction positive']), max(ref['fraction positive']))*1.1)
            # ax.legend()

            fig.savefig(os.path.join(self.working_dir, '%s_sample_%d.png' % (self.site_name, sample_index[n])))
            plt.close(fig)

        return compare_results


    @classmethod
    def plot_comparison(cls, fig, data, **kwargs):

        return None
