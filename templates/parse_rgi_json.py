#!/usr/bin/env python3

"""
Blablabla
"""

import os
import json 
import fnmatch
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import six


__version__ = "0.0.1"
__build__ = "03.12.2020"
__template__ = "PROCESS_RGI_FASTA-nf"

if __file__.endswith(".command.sh"):
    JSON_REPORTS = "$JSON_FILES".split()
    print("Running {} with parameters:".format(
        os.path.basename(__file__)))
    print("JSON_REPORTS: {}".format(JSON_REPORTS))


def render_mpl_table(data, col_width=3.0, row_height=0.625, font_size=14,
                     header_color='#40466e', row_colors=['#f1f1f2', 'w'], edge_color='w',
                     bbox=[0, 0, 1, 1], header_columns=0,
                     ax=None, **kwargs):
    if ax is None:
        size = (np.array(data.shape[::-1]) + np.array([0, 1])) * np.array([col_width, row_height])
        fig, ax = plt.subplots(figsize=size)
        ax.axis('off')

    mpl_table = ax.table(cellText=data.values, bbox=bbox, colLabels=data.columns, **kwargs)

    mpl_table.auto_set_font_size(False)
    mpl_table.set_fontsize(font_size)

    for k, cell in six.iteritems(mpl_table._cells):
        cell.set_edgecolor(edge_color)
        if k[0] == 0 or k[1] < header_columns:
            cell.set_text_props(weight='bold', color='w')
            cell.set_facecolor(header_color)
        else:
            cell.set_facecolor(row_colors[k[0]%len(row_colors) ])
    return ax


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
    

    # save df to image for report
    print(df_count_hits)

    ax = render_mpl_table(df_count_hits, header_columns=0, col_width=2.0)
    fig = ax.get_figure()
    fig.savefig('count_hits.png')

if __name__ == '__main__':
    main(JSON_REPORTS)