# -*- coding: utf-8 -*-
"""
Created on Wed Jan 03 10:48:31 2018

@author: jgarcia
"""

import matplotlib
matplotlib.use('Agg')

from matplotlib.collections import BrokenBarHCollection
import pandas as pd

### read ideogram:
file_name='example_chr1.pkl'
ideo = pd.read_pickle(file_name)

chromosome_list= ideo.chrom.unique()
# Height of each ideogram
chrom_height = .5

# Spacing between consecutive ideograms
chrom_spacing = .05

# Height of the gene track. Should be smaller than `chrom_spacing` in order to
# fit correctly
gene_height = 0.0

# Padding between the top of a gene track and its corresponding ideogram
gene_padding = 0.0


# Keep track of the y positions for ideograms and genes for each chromosome,
# and the center of each ideogram (which is where we'll put the ytick labels)
ybase = 0
chrom_ybase = {}
gene_ybase = {}
chrom_centers = {}

# Iterate in reverse so that items in the beginning of `chromosome_list` will
# appear at the top of the plot
for chrom in chromosome_list[::-1]:
    chrom_ybase[chrom] = ybase
    chrom_centers[chrom] = ybase + chrom_height / 2.
    gene_ybase[chrom] = ybase - gene_height - gene_padding
    ybase += chrom_height + chrom_spacing


layout = {'width': 2000, 'height': 1000, 'autosize': False, 'hovermode': 'closest', 'margin': {'l': 250, 'r': 199, 't': 120, 'b': 109, 'pad': 0}, 'xaxis1': {'anchor': 'y1', 'zeroline': False, 'ticks': 'inside', 'type': 'linear', 'range': [-2161775.8500000001, 45397292.850000001], 'showgrid': False, 'domain': [0.0, 1.0], 'side': 'bottom', 'tickfont': {'size': 10.0}, 'tick0': 0, 'dtick': 2000000, 'tickmode': False, 'mirror': 'ticks', 'showline': True}, 'yaxis1': {'anchor': 'x1', 'zeroline': False, 'ticks': 'inside', 'type': 'linear', 'range': [-1.0975000000000008, 23.047500000000014], 'showgrid': False, 'domain': [0.0, 1.0], 'side': 'left', 'tickfont': {'size': 10.0}, 'tick0': 21.700000000000014, 'dtick': -0.55000000000000071, 'tickmode': False, 'mirror': 'ticks', 'showline': True}, 'showlegend': False}

def return_figure():
    layout['shapes'] = []

    for chrom,group in ideo.groupby('chrom'):
        if chrom not in chromosome_list:
            continue
        for cramp in [x for x in range(600)]:
            layout['shapes'].append(
            {
            'type': 'rect',
            'x0': group.iloc[cramp,:].start,
            'x1': group.iloc[cramp,:].end,
            'y0': chrom_ybase[chrom],
            'y1': chrom_ybase[chrom] + chrom_height,
            'fillcolor': group.iloc[cramp,:].gieStain,
            'line': {
                'width': 0
            }
            }
            )
    example_figure = {
    'data': [],
    'layout': layout
    }
    return example_figure
#


#pp = pprint.PrettyPrinter(indent=4)
#pp.pprint(plotly_fig['layout'])



#plt.savefig('Ideo_' + Subject + '.png',bbox_inches = 'tight')
import flask
import dash
from dash.dependencies import Input, Output, State, Event
import dash_core_components as dcc
import dash_html_components as html

import pandas as pd

app = dash.Dash(__name__)
server = app.server

app.css.append_css({'external_url': 'https://codepen.io/chriddyp/pen/dZVMbK.css'})



app.layout = html.Div([
    
    html.Div([
    html.H3(
    id= "header1",className= "six columns",children= "Hello again"
    )],className= "row"),
    
    
#    html.Img(src='data:image/png;base64,{}'.format(encoded_image.decode()))
    html.Div([
    dcc.Graph(id= "ideogram",figure=return_figure())
    ]),
    
    ])

# Run the Dash app
if __name__ == '__main__':
    app.server.run(debug=True, threaded=True)
