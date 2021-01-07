#!/usr/bin/env python3

import ast
import pandas as pd
import plotly.graph_objs as go
from plotly.offline import plot
from plotly.subplots import make_subplots

SUMMARY_TABLE = "$summary_df"
HITS_TABLE = "$hit_table"
HEATMAP = "$heatmap_df"

colorscale = ['#f0f0f0', '#5a6e92','#04a0fc']

conversion = {0: "No hit", 1: "Strict", 2: "Perfect"}

def convert_tuple(tup):
    try:
        tup = ast.literal_eval(tup)
        str =  ', '.join(tup) 
        return str
    except:
        return tup

def df_to_plotly(df):
    return {'z': df.values.tolist(),
            'x': df.columns.tolist(),
            'y': df.index.tolist()}

def get_hovertext(df):
    hovertext = list()
    for yi, yy in enumerate(df.index):
        hovertext.append(list())
        for xi, xx in enumerate(df.columns.tolist()):
            hovertext[-1].append('Sample: {}<br />Symbol: {}<br />Hit Type: {}'.format(xx, yy, conversion[df.iat[yi, xi]]))
    return hovertext
    
def main(summary_file, hits_table, heatmap):

    # load data files
    #   - summary: number of hits per sample 
    summary_df = pd.read_csv(summary_file)
    #   - heatmap: pivot table data
    df_heatmat = pd.read_csv(heatmap)
    df_heatmat = df_heatmat.set_index(df_heatmat.columns[0])
    #   - hits table: comprehensive results table
    hits_df = pd.read_csv(hits_table)
    hits_df = hits_df.applymap(convert_tuple)

    # main report skeleton 
    fig = make_subplots(rows=11, cols=2, vertical_spacing=0.1, 
                    specs=[[{"type": "table", "rowspan": 1, "colspan": 2}, None],
                           [None, None],
                           [{"type": "table", "rowspan": 2, "colspan": 2}, None],
                           [None, None],
                           [{"type": "xy", "rowspan":4, "colspan": 2}, None],
                           [None, None],
                           [None, None],
                           [None, None],
                           [{"type": "table", "rowspan":3 ,"colspan": 2}, None],
                           [None, None],
                           [None, None]],
                        subplot_titles=("","Summary Table", "Hits Heatmap", "Hits Table"))

    # intro text
    summary = ('RGI-nf integrates <a href="https://github.com/arpcard/rgi">The Resistance Gene Identifier (RGI)</a> package<br>'
               'to predict resistome(s) from nucleotide data based on homology and SNP models.<br> '
               ' <br>'
               'The application uses reference data from the <a href="https://card.mcmaster.ca/">Comprehensive Antibiotic Resistance Database (CARD)</a>.'
               ' <br>'
               'A total of {0} samples were analysed, detecting a total of {1} unique AMR genes.'
               ' <br>'.format(len(summary_df.columns[0]), len(hits_df.index)))

    fig.add_annotation(x=0, xref='paper', y=1, yref='paper', text=summary,
                       showarrow=False, font=dict(size=14), align='center', 
                       bgcolor='#f0f0f0')

    # summary table
    fig.add_trace(go.Table(header=dict(values=["Sample", "Perfect", "Strict", "Loose"],
                           font=dict(size=14), align="center", fill_color="#04a0fc"),
                           cells=dict(values=[summary_df[k].tolist() for k in summary_df.columns[0:]], 
                                      align = ["left", "center"], fill_color='#f0f0f0')),
                  row=3, col=1)
    
    # heatmap 
    hovertext = get_hovertext(df_heatmat)
    fig.add_trace(go.Heatmap(df_to_plotly(df_heatmat), colorscale=colorscale, 
                             hoverinfo='text', text=hovertext, showscale=False), row=5, col=1)
    
    # hits table
    fig.add_trace(go.Table(header=dict(values=["Gene Symbol", "Gene Family", "Drug Class", "Resistance Mechanism"],
                        font=dict(size=14), align="center", fill_color="#04a0fc"),
                        cells=dict(values=[hits_df[k].tolist() for k in hits_df.columns[0:]], 
                                    align = "left",  fill_color='#f0f0f0')),
                row=9, col=1)
    
    # title
    fig.update_layout(title={'text': "<b>RGI-nf Report</b>",
                            'y':0.95,
                            'x':0.5,
                            'xanchor': 'center',
                            'yanchor': 'top',
                            'font': {'size': 24, 'color': '#04a0fc'}})
    fig.update_layout(height=1080, template='ggplot2', plot_bgcolor='rgba(0,0,0,0)')
    fig.write_html("multiqc_report.html")
    

if __name__ == "__main__":
    main(SUMMARY_TABLE, HITS_TABLE, HEATMAP)
