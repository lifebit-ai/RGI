#!/usr/bin/env python3

import io
import os
import json 
import fnmatch
import pandas as pd
import numpy as np
from pretty_html_table import build_table
import matplotlib
from collections import defaultdict, Counter
import seaborn as sns
import matplotlib.pyplot as plt

__version__ = "0.0.1"
__build__ = "17.12.2020"
__template__ = "PROCESS_RGI_BWT-nf"

if __file__.endswith(".command.sh"):
    JSON_REPORTS = "$JSON_FILES".split()
    print("Running {} with parameters:".format(
        os.path.basename(__file__)))
    print("JSON_REPORTS: {}".format(JSON_REPORTS))

def main(json_reports):

    # main hit counts df 
    df_hits = pd.DataFrame(columns=['Sample', 'Gene Symbol', 'Hit Type', 'Gene Family', 'Drug Class', 'Resistance Mechanism'])

    df_count_hits = pd.DataFrame(columns=['Sample', 'Perfect', 'Strict', 'Loose'])

    for output_json_file in json_reports:

        print("Processing count hits file: {}".format(output_json_file))
        sample_name = os.path.basename(output_json_file).split('_')[0]
        perfect = 0
        sctrict = 0
        loose = 0

        with open(output_json_file) as json_fh:
            report_json = json.load(json_fh)
            for hit in report_json:
                hit = dict(hit)
                gene_symbol = hit['cvterm_name']
                if int(hit['length_coverage']['uncovered']) == 0:
                    hit_type = 'Perfect'
                    perfect += 1
                elif int(hit['length_coverage']['covered']) >= (int(hit['reference']['sequence_length']) * 0.8):
                     hit_type = 'Strict'
                     sctrict += 1
                else:
                    hit_type = 'Loose'
                    loose += 1
                gene_family = hit['resistomes']['AMR Gene Family']
                drug_class = hit['resistomes']['Drug Class']
                resistance_mechanism = hit['resistomes']['Resistance Mechanism']

                df_hits = df_hits.append({
                            'Sample': sample_name,
                            'Gene Symbol': gene_symbol, 
                            'Hit Type': hit_type, 
                            'Gene Family': gene_family,
                            'Drug Class': drug_class,
                            'Resistance Mechanism': resistance_mechanism}, ignore_index=True)
        df_count_hits = df_count_hits.append({
            'Sample': sample_name,
            'Perfect': perfect,
            'Strict': sctrict,
            'Loose': loose
        }, ignore_index=True)
        

    df_count_hits.to_csv("results_summary.csv", index=False)

    """
    html_str = build_table(df_count_hits, "grey_dark")
    
    with open("results_summary.html", "w") as html_fh:
        html_fh.write(html_str)
    """
    
    # pivot dataframe
    genes_df = pd.pivot(df_hits, columns='Sample', index='Gene Symbol', values='Hit Type')
    #print(genes_df)

    # Create a dictionary that will convert type of hit to num. value
    genes_df = genes_df.replace({"Perfect": 2, "Strict": 1, "Loose":0})

    # drop rows with all zeros in all columns
    genes_df = genes_df.loc[~(genes_df==0).all(axis=1)]
    print(genes_df)
    genes_df.to_csv("card_hits.csv")

    # Fixed colourmap values (white, gray-blue, bright blue)
    cmap_values = [0, 1, 2, 3]
    custom_cmap = matplotlib.colors.ListedColormap(['#ffffff', '#5a6e92', '#04a0fc'])
    norm = matplotlib.colors.BoundaryNorm(cmap_values, custom_cmap.N)

    # create plot
    sns.heatmap(genes_df, cmap=custom_cmap, cbar=False, norm=norm)
    sns.set_style("white")
    plt.savefig("card_hits_headmap.png", bbox_inches="tight", format="png", pad_inches=0.5)


    # create table with classification
    results_table = df_hits.drop(columns=['Sample', 'Hit Type'])
    results_table.to_csv("results_hits.csv", index=False)

    """
    html_str = build_table(results_table, "grey_dark")
    with open("results_hits.html", "w") as html_fh:
        html_fh.write(html_str)
    """

if __name__ == '__main__':
    main(JSON_REPORTS)
    #main(["/Users/inesmendes/lifebit/git/RGI/work/62/648a6c202c8e490e96f6b7fbe793b4/ERR4008003_rgi_bwt.allele_mapping_data.json"])