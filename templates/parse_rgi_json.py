#!/usr/bin/env python3

"""
Blablabla
"""

import io
import os
import json 
import fnmatch
import pandas as pd
import numpy as np
from pretty_html_table import build_table

__version__ = "0.0.1"
__build__ = "03.12.2020"
__template__ = "PROCESS_RGI_FASTA-nf"

if __file__.endswith(".command.sh"):
    JSON_REPORTS = "$JSON_FILES".split()
    print("Running {} with parameters:".format(
        os.path.basename(__file__)))
    print("JSON_REPORTS: {}".format(JSON_REPORTS))


def main(json_reports):

    # main hit counts df 
    df_count_hits = pd.DataFrame(columns=['Sample', 'Perfect', 'Strict', 'Loose'])

    for count_hits_json_file in json_reports:
        print("Processing count hits file: {}".format(count_hits_json_file))
        with open(count_hits_json_file) as counts_fh:
            count_hits_json = json.load(counts_fh)
            df_count_hits = df_count_hits.append({
                            'Sample': count_hits_json_file.replace("_card_rgi_parsed-count-hits.json", ''),
                            'Perfect': int(count_hits_json['Perfect']), 
                            'Strict': int(count_hits_json['Strict']), 
                            'Loose': int(count_hits_json['Loose'])}, ignore_index=True)
    
    print(df_count_hits)

    # save df to html table for report
    #str_io = io.StringIO()
    #df_count_hits.to_html(buf=str_io, classes='table table-striped')
    #html_str = str_io.getvalue()
    html_str = build_table(df_count_hits, "grey_dark")
    
    with open("results_summary.html", "w") as html_fh:
        html_fh.write(html_str)

if __name__ == '__main__':
    main(JSON_REPORTS)