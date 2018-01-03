# -*- coding: utf-8 -*-
"""
Created on Wed Jan 03 10:48:31 2018

@author: jgarcia
"""

import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt
from matplotlib.collections import BrokenBarHCollection
import pandas as pd

import pprint
import numpy as np

# Plotly
import plotly.plotly as py
import plotly.tools as tls

def chromosome_collections(df, y_positions, height,  **kwargs):
    """
    Yields BrokenBarHCollection of features that can be added to an Axes
    object.
    Parameters
    ----------
    df : pandas.DataFrame
        Must at least have columns ['chrom', 'start', 'end', 'color']. If no
        column 'width', it will be calculated from start/end.
    y_positions : dict
        Keys are chromosomes, values are y-value at which to anchor the
        BrokenBarHCollection
    height : float
        Height of each BrokenBarHCollection
    Additional kwargs are passed to BrokenBarHCollection
    """
    del_width = False
    if 'width' not in df.columns:
        del_width = True
        df['width'] = df['end'] - df['start']
    for chrom, group in df.groupby('chrom'):
        print(chrom)
        yrange = (y_positions[chrom], height)
        xranges = group[['start', 'width']].values
        yield BrokenBarHCollection(
            xranges, yrange, facecolors=group['colors'], **kwargs)
    if del_width:
        del df['width']


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

###
# Colors for different chromosome stains
color_lookup = {
    'red': [255, 0, 0],
    'yellow': [255, 255, 0],
    'blue': [0, 0, 255],
    'orange': [255, 165, 0],
    'green': [50, 205, 50],
    'black': [0, 0, 0],
    'purple': [128, 0, 128],
    'silver': [211, 211, 211],
}


ideo['colors'] = ideo['gieStain'].apply(lambda x: tuple([round(y / float(255),1) for y in color_lookup[x]]))
# Add a new column for width
ideo['width'] = ideo.end - ideo.start

# Width, height (in inches)
figsize = (20, 10)

fig = plt.figure(figsize=figsize)
ax = fig.add_subplot(111)

# Now all we have to do is call our function for the ideogram data...
print("adding ideograms...")
for collection in chromosome_collections(ideo, chrom_ybase, chrom_height, edgecolors=None, linewidths= 0):
    ax.add_collection(collection)


# Axes tweaking
ax.set_xticks([x for x in range(0,max(ideo.end),2000000)])
plt.xticks(fontsize = 10,rotation = 90)
ax.tick_params(axis = 'x',pad = 10)

ax.tick_params(axis='y', which='major', pad=30)
ax.set_yticks([chrom_centers[i] for i in chromosome_list])
ax.set_yticklabels(chromosome_list, fontsize = 10)
ax.axis('tight')

from mpld3 import fig_to_dict

test= fig_to_dict(fig)
#print(test.keys())
#print(test['axes'][0])
#print(test['id'])
#print (test)


#from mplexporter.utils import color_to_hex
#from mplexporter.exporter import Exporter
#from mplexporter.renderers import Renderer
#

#renderer = PlotlyRenderer()
#exporter = Exporter(renderer)
#exporter.run(fig)
#print(renderer.data)

plotly_fig= tls.mpl_to_plotly(fig, strip_style=False, verbose=False)
print(tls.mpl_to_plotly(fig))

plotly_fig['layout']['shapes'] = []

for chrom,group in ideo.groupby('chrom'):
    for cramp in [x for x in range(500,800)]:
        plotly_fig['layout']['shapes'].append(
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

#

#pp = pprint.PrettyPrinter(indent=4)
#pp.pprint(plotly_fig['layout'])


#plt.savefig('Ideo_' + Subject + '.png',bbox_inches = 'tight')
import flask
import dash
from dash.dependencies import Input, Output, State, Event
import dash_core_components as dcc
import dash_html_components as html
import base64


image_filename = 'Ideo_IRIS_313-11825.png' # replace with your own image
encoded_image = base64.b64encode(open(image_filename, 'rb').read())


import pandas as pd

app = dash.Dash(__name__)
server = app.server

app.css.append_css({'external_url': 'https://codepen.io/chriddyp/pen/dZVMbK.css'})



app.layout = html.Div([
    
    
#    html.Img(src='data:image/png;base64,{}'.format(encoded_image.decode()))
    html.Div([
    dcc.Graph(id= "ideogram",figure=plotly_fig)
    ]),
    
    ])

# Run the Dash app
if __name__ == '__main__':
    app.server.run(debug=True, threaded=True)
