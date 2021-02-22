import matplotlib.pyplot as plt
from matplotlib import cm
import seaborn as sns
import sqlite3
import sys, os, json, collections, struct, datetime
import pandas as pd
import numpy as np
import matplotlib as mpl
import os


if __name__ == '__main__':

    exp_name = r'FPG_longrun_MAF_119_6_FIXED_MSP_5years_pop_scale_1.0_v2'
    wdir = fr'C:\Users\jorussell\Dropbox (IDM)\Malaria Team Folder\projects\parasite_genetics\DTK\longrun_genetic_behavior\{exp_name}'
    #iterate through all sims in this dir
    for dirpath,dirname,filenames in os.walk(wdir):
        if len(filenames) == 0:
            print(f'CAUTION: no_reports_in_dir_level_{dirpath}')
            continue
        else:
            sql_report_fn = os.path.join(dirpath,filenames[1])
            conn = sqlite3.connect(sql_report_fn)
            cursor = conn.cursor()
            sql = ("SELECT ID.InfectionID,ID.SimTime, INF.HumanID, INF.GenomeID, GSD.AlleleRoots, GSD.NucleotideSequence, GSD.GenomeLocation FROM"
                   "(SELECT InfectionID, SimTime "
                   "FROM InfectionData ) ID "
                   "JOIN "
                   "(SELECT InfectionID, HumanID, GenomeID "
                   "FROM Infections) INF "
                   "ON (ID.InfectionID = INF.InfectionID) "
                   "JOIN "
                   "(SELECT GenomeID, AlleleRoots, NucleotideSequence, GenomeLocation "
                   "FROM GenomeSequenceData ) GSD "
                   "ON (INF.GenomeID = GSD.GenomeID) "
                   "WHERE (ID.SimTime BETWEEN 100 AND 105) AND (GSD.GenomeLocation != 200000) "
                   )
            cursor.execute(sql)
# within day what is pairwise IBD of sampled genomes as SELECT
#
            res = cursor.fetchall()
            df = pd.DataFrame(res, columns=[x[0] for x in cursor.description])
            df2 = df.groupby(['SimTime', 'HumanID', 'InfectionID'])['AlleleRoots'].apply(list).reset_index(name='AlleleRoots')
            df3 = df.groupby(['SimTime', 'HumanID', 'InfectionID'])['NucleotideSequence'].apply(list).reset_index(
                name='NucleotideSequence')

            df3['NucleotideSequence'] = pd.Series([[1 if elem == 3 else 0 for elem in l] for l in df3['NucleotideSequence']])


            COI = df3.groupby(['SimTime','HumanID'])['InfectionID'].apply(list)
            # save df to file in dirpath

            conn.close()
    print('success!')