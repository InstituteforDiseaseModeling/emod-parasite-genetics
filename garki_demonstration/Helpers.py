import pandas as pd
import json


def reference_prevalence(metadata):

    with open(metadata['fname']) as json_file:
        data = json.load(json_file)

        true_pfpr = data['Channels']['True Prevalence']['Data'][metadata['start_date']:(metadata['start_date']+metadata['analyze_duration'])]

        df = pd.DataFrame({'True_PfPR': true_pfpr})
        # df = df.set_index('True_PfPR')

    return df