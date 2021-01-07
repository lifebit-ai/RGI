#!/usr/bin/env python3

import pandas as pd
import plotly.graph_objs as go
from plotly.offline import plot
from plotly.subplots import make_subplots

SUMMARY_TABLE = "$summary_df"
HITS_TABLE = "$hit_table"
HEATMAP = "$heatmap_df"

def create_table_tracer(header_values, header_font, header_line, header_fill,
                        cells_values, cells_font, cells_line, cells_fill,
                        domain):
    """
    """

    tracer = go.Table(header=dict(values=header_values,
                                  font=header_font,
                                  align='left',
                                  line=header_line,
                                  fill=header_fill),
                      cells=dict(values=cells_values,
                                 font=cells_font,
                                 align='left',
                                 line=cells_line,
                                 fill=cells_fill),
                      domain=domain)

    return tracer


def df_to_plotly(df):
    return {'z': df.values.tolist(),
            'x': df.columns.tolist(),
            'y': df.index.tolist()}


def main(summary_file, hits_table, heatmap):
    
    # main report 
    fig = make_subplots(rows=3, cols=2, vertical_spacing=0.1, 
                        specs=[[{"type": "table"}, None],
                              [{"colspan": 2}, None],
                              [{"colspan": 2, "type": "table"}, None]],
                        subplot_titles=("Summary Table", "Hits Heatmap", "Hits Table"))

    # summary table
    summary_df = pd.read_csv(summary_file)
    fig.add_trace(go.Table(header=dict(values=["Sample", "Perfect", "Strict", "Loose"],
                           font=dict(size=10), align="left"),
                           cells=dict(values=[summary_df[k].tolist() for k in summary_df.columns[0:]], 
                                      align = "left")),
                  row=1, col=1)
    
    # heatmap 
    conversion = {0: "No hit", 1: "Strict", 2: "Perfect"}

    df_heatmat = pd.read_csv(heatmap)
    df_heatmat = df_heatmat.set_index(df_heatmat.columns[0])
    colorscale = ['#ffffff', '#5a6e92','#04a0fc']
    hovertext = list()
    for yi, yy in enumerate(df_heatmat.index):
        hovertext.append(list())
        for xi, xx in enumerate(df_heatmat.columns.tolist()):
            hovertext[-1].append('Sample: {}<br />Symbol: {}<br />Hit Type: {}'.format(xx, yy, conversion[df_heatmat.iat[yi, xi]]))
    
    fig.add_trace(go.Heatmap(df_to_plotly(df_heatmat), colorscale=colorscale, 
    hoverinfo='text', text=hovertext, showscale=False), row=2, col=1)

    fig.update_layout(title={'text': "RGI-nf Report",
                             'y':0.95,
                             'x':0.5,
                             'xanchor': 'center',
                             'yanchor': 'top',
                             'font': {'size': 24}})
    
    # hits table
    hits_df = pd.read_csv(hits_table)
    fig.add_trace(go.Table(header=dict(values=["Gene Symbol", "Gene Family", "Drug Class", "Resistance Mechanism"],
                        font=dict(size=10), align="left"),
                        cells=dict(values=[hits_df[k].tolist() for k in hits_df.columns[0:]], 
                                    align = "left")),
                row=3, col=1)
    fig.show()
    

if __name__ == "__main__":
    #main(SUMMARY_TABLE, HITS_TABLE, HEATMAP)
    main("../work/d2/d00c5520af3d7ee4144cdb6cb4bf61/results_summary.csv", 
         "../work/d2/d00c5520af3d7ee4144cdb6cb4bf61/results_hits.csv", 
         "../work/9d/61c13da64d11bcd7ede06d8941e143/card_hits.csv")