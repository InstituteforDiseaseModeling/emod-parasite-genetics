import json
import os
from io import BytesIO

from simtools.Analysis.AnalyzeManager import AnalyzeManager
from simtools.Analysis.BaseAnalyzers import BaseAnalyzer


class DownloadAnalyzer_batch(BaseAnalyzer):
    """
    This analyzer is based on the DownloadAnalyzer and allows the download of files based on tags
    """

    def __init__(self, filenames, output_path="output"):
        super().__init__(filenames=filenames, parse=False)

        self.output_path = output_path
        self.filenames = filenames
    def initialize(self):
        super().initialize()
        if not os.path.exists(self.output_path):
            os.makedirs(self.output_path,exist_ok=True)

    def per_experiment(self, experiment):
        """
        Set and create the output path.
        :param experiment: experiment object to make output directory for
        :return: Nothing
        """
        exp_name = experiment.exp_name
        output_dir = fr'C:\Users\jorussell\Dropbox (IDM)\Malaria Team Folder\projects\parasite_genetics\DTK\longrun_genetic_behavior\{exp_name}_v2'
        self.output_path = output_dir

        if not os.path.exists(self.output_path):
            os.makedirs(self.output_path,exist_ok=True)


    def get_sim_folder(self, sim_id):
        """
        Concatenate the specified top-level output folder with the simulation ID
        :param parser: A simulation output parsing thread
        :return: The name of the folder to download this simulation's output to
        """
        return os.path.join(self.output_path, sim_id)

    def select_simulation_data(self, data, simulation):
        # Create a folder for the current simulation
        lhm_tag = str(simulation.tags['larval_habitat_multiplier'])[2:5]
        if lhm_tag == 0:
            lhm_tag =1
        seed = simulation.tags['Run_Number']
        sim_folder = os.path.join(self.output_path,f"LHM_{lhm_tag}_seed_{seed}")
        # sim_folder = os.path.join(self.output_path,f"seed_{seed}")

        if not os.path.exists(sim_folder):
            os.makedirs(sim_folder, exist_ok=True)
        # Create the requested files
        for filename in self.filenames:
            file_path = os.path.join(sim_folder, os.path.basename(filename))

            with open(file_path, 'wb') as outfile:
                outfile.write(data[filename])


if __name__ == '__main__':

    REPORT_NAME = 'InsetChart.json'
    # REPORT_NAME = 'MalariaSqlReport.db'

    exp_ids = ['c0c8aaef-3f6f-eb11-a2dd-c4346bcb7271']
    # exp_ids = ['194e33f2-936d-eb11-a2dd-c4346bcb7271']

    analyzer = DownloadAnalyzer_batch(filenames=[f'output\\{REPORT_NAME}'])


    am = AnalyzeManager(exp_list=exp_ids, analyzers=analyzer)
    am.analyze()
