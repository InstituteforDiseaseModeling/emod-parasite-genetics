
import matplotlib.pyplot as plt
from matplotlib import cm
import seaborn as sns
import sqlite3

import sys, os, json, collections, struct, datetime
import pandas as pd
import numpy as np

import matplotlib as mpl
from simtools.Analysis.AnalyzeManager import AnalyzeManager
from simtools.Analysis.BaseAnalyzers import BaseAnalyzer
from simtools.SetupParser import SetupParser
from scipy import interpolate
import os

wdir = os.getcwd()

def CompareValues(messages, var, exp, act):
    """
    Compare two values - My unit test like feature
    """
    success = True
    if (abs(exp - act) > 0.005):
        messages.append(var + ": Expected " + str(exp) + " but got " + str(act))
        success = False
    return success


def CompareArrays(messages, name, exp_data, act_data):
    """
    Compare to the two arrays to make sure the values are 'very' similar
    """

    messages.append("!!!!!!!!!!!!!!!!!! Compare " + name + " !!!!!!!!!!!!!!!!!!!!!!!")
    success = True

    exp_num = len(exp_data)
    act_num = len(act_data)
    success = CompareValues(messages, name + ":num", exp_num, act_num)
    if not success:
        return success

    i = 0
    for exp_val, act_val in zip(exp_data, act_data):
        # print(str(i)+"-"+str(exp_val)+" ? "+str(act_val))
        success = CompareValues(messages, name + ":val[" + str(i) + "]", exp_val, act_val)
        if not success:
            return success
        i += 1

    messages.append("!!!!!!!!!!!!!!!!!! PASSED - " + name + " !!!!!!!!!!!!!!!!!!!!!!!")

    return success


class Campaign:
    """
    A class for reading an campaign.json file.
    """

    def __init__(self):
        self.fn = ""
        self.json_data = {}

    def Read(self, filename, messages):
        """
        Read the file given by filename and verify that the file has two CampaignEvents.
        """
        self.fn = filename

        with open(filename, 'r') as file:
            self.json_data = json.load(file)

        exp_num_events = 2  # one for outbreak and one for drugs
        act_num_events = len(self.json_data["Events"])
        CompareValues(messages, "Campaign:Num Events", exp_num_events, act_num_events)

    def GetBarcodeAlleleFrequencies(self, messages):
        """
        The outbreak is expected to be in the first event.
        """
        freqs = self.json_data["Events"][0]["Event_Coordinator_Config"]["Intervention_Config"][
            "Barcode_Allele_Frequencies_Per_Genome_Location"]
        print(freqs)
        return freqs


class MalariaSqlReport:
    """
    A class for reading MalariaSqlReport
    """

    def __init__(self):
        self.fn = ""
        self.conn = None
        self.cursor = None

    def Open(self, filename):
        """
        Open the database specified by filename and prepare for queries
        """
        self.fn = filename
        self.conn = sqlite3.connect(self.fn)
        self.conn.row_factory = lambda cursor, row: row[0]  # makes the output lists of values instead of tuples
        self.cursor = self.conn.cursor()

    def Close(self):
        self.conn.close()

    def GetDataAsInsetChartChannel(self, channel_name):
        """
        Extract data from the database and format it into a list of values similar to what
        is in the corresponding InsetChart channel.
        """
        if channel_name == "Statistical Population":
            self.cursor.execute("SELECT COUNT(*) FROM Health GROUP BY SimTime")
            return self.cursor.fetchall()
        elif channel_name == "Infected":
            self.cursor.execute("SELECT COUNT(*) FROM Health GROUP BY SimTime")
            pops = self.cursor.fetchall()
            self.cursor.execute("SELECT SUM(Infected) " \
                                "FROM (SELECT Health.SimTime, 0 as Infected " \
                                "FROM Health " \
                                "GROUP BY Health.SimTime " \
                                "UNION " \
                                "SELECT SimTime, COUNT(*) AS Infected " \
                                "FROM (SELECT InfectionData.SimTime AS SimTime, Infections.HumanID AS HumanID " \
                                "FROM Infections INNER JOIN InfectionData ON InfectionData.InfectionID = Infections.InfectionID " \
                                "GROUP BY SimTime, HumanID " \
                                "ORDER BY SimTime) " \
                                "GROUP BY SimTime) " \
                                "GROUP BY SimTime")
            infs = self.cursor.fetchall()
            infected = [int(i) / int(p) for i, p in zip(infs, pops)]
            return infected

def application(output_path="output"):
    print("!!!!! Check allele frequencies !!!!!")

    # Define the names of the files to be used in the test
    campaign_fn = "campaign.json"
    sql_db_fn = os.path.join(output_path, "MalariaSqlReport.db")

    messages = []

    # Open all of the files
    campaign = Campaign()
    campaign.Read(campaign_fn, messages)

    sql_db = MalariaSqlReport()
    sql_db.Open(sql_db_fn)

    # Get the expected allele frequencies from the campaign file
    exp_allele_freqs = campaign.GetBarcodeAlleleFrequencies(messages)

    # Get the actual allele frequencies from the database
    act_allele_freqs = sql_db.GetBarcodeAlleleFrequencies(messages)

    # Verify the reports have similar data
    CompareValues(messages, "Campaign vs SqlDB:Num Allele Frequencies", len(exp_allele_freqs), len(act_allele_freqs))
    for i in range(len(exp_allele_freqs)):
        CompareArrays(messages, "Campaign vs SqlDB:Allele Frequencies[" + str(i) + "]", exp_allele_freqs[i],
                      act_allele_freqs[i])

    sql_db.Close()

    # Create a file with the results of the test
    output = {}
    output["messages"] = []
    for line in messages:
        output["messages"].append(line)

    report_fn = os.path.join(output_path, "report_sync_check.json")
    with open(report_fn, 'w') as report_file:
        json.dump(output, report_file, indent=4)

    for line in messages:
        print(line)

    print("!!!!! Done checking allele frequencies !!!!!")


class Genetics_Analyzer(BaseAnalyzer):
    def __init__(self, output_fname):
        super(Genetics_Analyzer, self).__init__()
        self.filenames = ['output/InsetChart.json','output/MalariaSqlReport.db']
        self.output_fname = output_fname


    def select_simulation_data(self, data, simulation):

        simdata = pd.DataFrame()
        sql_db = MalariaSqlReport()
        sql_db.Open(filename = r'C:\Users\jorussell\Downloads\MalariaSqlReport_test_62.db')
        sql_db.GetDataAsInsetChartChannel(channel_name='Infected')

        sql_db.Close()

        return simdata

    def finalize(self, all_data):

        selected = [data for sim, data in all_data.items()]
        if len(selected) == 0:
            print("No data have been returned... Exiting...")
            return



if __name__ == '__main__' :


    expids = ['ad854768-8c65-eb11-a2dd-c4346bcb7271']


    expnames = ['base_demo']
    channel_name = 'PCR Parasite Prevalence'
    for expname, expid in zip(expnames, expids) :
        output_fname = os.path.join(wdir, expid,expname)
        am = AnalyzeManager(exp_list = expids,
                            analyzers=Genetics_Analyzer(output_fname))
        am.analyze()