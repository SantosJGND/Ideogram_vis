# -*- coding: utf-8 -*-
"""
Created on Fri Jan 05 16:33:18 2018

@author: jgarcia
"""


### imports
import pandas as pd
import dash
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objs as go
import plotly.figure_factory as ff
import sqlite3
import base64
import numpy as np
#import pandas as pd
import collections
#import copy
#
from sklearn import preprocessing
from sklearn.decomposition import PCA


def recursively_default_dict():
        return collections.defaultdict(recursively_default_dict)


###

def db_profile_query(filename,Region):
    ret= []
    try:
        con= sqlite3.connect(filename)
        cur= con.cursor()
        
        query= 'SELECT Chr, window, cluster, profiles, blocks FROM Profiles WHERE window >= ? AND window < ?'
#        print(query.format(Region))
        cur.execute(query,[str(x) for x in Region])
#        cur.execute('SELECT * FROM Profiles;')
        data= cur.fetchall()
        for line in data:
            ret.append([line[0],line[1],[float(x) for x in line[3].split(';')],[int(x) for x in line[4].split(';')]])
        
        cur.execute('SELECT blocks FROM Profiles WHERE window = -1')
        data= cur.fetchone()
        Names= data[0].split(';')
        return ret, Names
        
    except sqlite3.Error:
        if con:
            con.rollback()
            print('something')
    finally:
        if con:
            con.close()



def db_blocks_query(filename,Region):
    ret= recursively_default_dict()
    try:
        con= sqlite3.connect(filename)
        cur= con.cursor()
        
        query= 'SELECT Chr, window, profile FROM Blocks WHERE window >= ? AND window < ?'
#        print(query.format(Region))
        cur.execute(query,[str(x) for x in Region])
#        cur.execute('SELECT * FROM Profiles;')
        data= cur.fetchall()
        for line in data:
            ret[line[0]][line[1]] = [int(c) for c in line[2].split(';')]
        
        cur.execute('SELECT profile FROM Blocks WHERE window = -1')
        data= cur.fetchone()
        Names= data[0].split(';')
        return ret, Names
    except sqlite3.Error:
        if con:
            con.rollback()
            print('something')
    finally:
        if con:
            con.close()


def generate_table(dataframe):
    return html.Table(
        # Header
        [html.Tr([html.Th(col) for col in dataframe.columns])] +

        # Body
        [html.Tr([
            html.Td(dataframe.iloc[i][col]) for col in dataframe.columns
        ]) for i in range(len(dataframe))]
    )


### function to return figures.
def return_figure(ideo,layout):
    
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
    
    
    layout['shapes'] = []

    for chrom,group in ideo.groupby('chrom'):
        if chrom not in chromosome_list:
            continue
        for cramp in [x for x in range(group.shape[0])]:
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


def loadings_graph(df,orderCore,scheme,coords,pop_refs,color_here,opac):
    figure = {
        'data': [go.Scatter3d(
            x = df.iloc[coords[i],df.columns.get_loc('PC1')],
            y = df.iloc[coords[i],df.columns.get_loc('PC2')],
            z = df.iloc[coords[i],df.columns.get_loc('PC3')],
            type='scatter3d',
            mode= "markers",
            text= orderCore.iloc[coords[i],:][["ID","NAME","COUNTRY","Initial_subpop"]].apply(lambda lbgf: (
          "<b>{}</b><br>Name: {}<br>Country: {}<br>{}".format(lbgf[0],lbgf[1],lbgf[2],lbgf[3])),
            axis= 1),
        marker= {
    #        'color': scheme,
            'color': color_here[i],
            'line': {'width': 0},
            'size': 4,
            'symbol': 'circle',
          "opacity": opac
          },
          name= pop_refs[i]
        ) for i in list(set(scheme))],
        'layout': {
      "autosize": True, 
      "hovermode": "closest", 
      "legend": {
      "x": 0.798366013072, 
      "y": 0.786064620514, 
      "borderwidth": 1, 
      "font": {"size": 13}
          },
      "margin": {
        "r": 40, 
        "t": 50, 
        "b": 20, 
        "l": 30, 
        "pad": 0
      }, 
      "scene": {
        "aspectmode": "auto", 
        "aspectratio": {
          "x": 1.02391505715, 
          "y": 0.737436541286, 
          "z": 1.3243763495
        }, 
        "camera": {
          "center": {
            "x": 0, 
            "y": 0, 
            "z": 0
          }, 
          "eye": {
            "x": 1.80578427889, 
            "y": 1.17729688569, 
            "z": 0.201532084509
          }, 
          "up": {
            "x": 0, 
            "y": 0, 
            "z": 1
          }
        }, 
        "xaxis": {
          "title": "PC1", 
          "type": "linear"
        }, 
        "yaxis": {
          "title": "PC2", 
          "type": "linear"
        }, 
        "zaxis": {
          "title": "PC3", 
          "type": "linear"
        }
      }, 
      "showlegend": True, 
      "title": "<b>Accessions - loadings</b>", 
      "xaxis": {"title": "V3"}, 
      "yaxis": {"title": "V2"}
    }
    }
    
    return figure
###
### corss reference Blocks and Profiles and Analyse to prodduce dataframes

def Clover(Profiles,focus_indexes,Trend,chromosomes,threshold,X_threshold,Focus,target,Region):
#        print(len(Profiles[1]))
#        print([x for x in Profiles[1].keys()][:10])
#        print('hey')
        #### Chose_profiles: automatically chose clusters with at least one included 
        #### focus accession of 'target' color at a threshold >= target_threshold
#        Chose_profiles = {CHR:{bl:[y for y in Profiles[CHR][bl].keys() if sum([int(Profiles[CHR][bl][y][z] >= X_threshold) \
#        for z in [x for x in focus_indexes if Blocks[CHR][bl][x] in target]]) >= 1] \
#        for bl in Profiles[CHR].keys() if bl >= min(Region)*1e6 and bl <= max(Region)*1e6 and\
#        len([x for x in focus_indexes if Blocks[CHR][bl][x] in target]) / float(len(Focus)) >= threshold} \
#        for CHR in Blocks.keys() if CHR in chromosomes}
        
        Clover= [cl[2] for cl in Profiles if len([x for x in range(len(cl[2])) if cl[2][x] >= X_threshold and x in focus_indexes and cl[3][x] in target]) > 1]
#    #    Coordinates = [[[[CHR,bl] for x in Chose_profiles[CHR][bl]] for bl in sorted(Chose_profiles[CHR].keys())] for CHR in sorted(Chose_profiles.keys())]
#    #    Coordinates = [z for z in it.chain(*[y for y in it.chain([x for x in it.chain(*Coordinates)])])]
#    #    Coordinates= np.array(Coordinates)
#        
#        Clover= [[[Profiles[CHR][bl][x] for x in Chose_profiles[CHR][bl]] for bl in sorted(Chose_profiles[CHR].keys())] for CHR in sorted(Chose_profiles.keys())]
#    #    del Blocks
#    #    del Profiles
#        Clover= [z for z in it.chain(*[y for y in it.chain(*Clover)])]
        Clover= np.array(Clover)
        Clover = np.nan_to_num(Clover)
        Test = Clover
        
        print('Clover shape: ', Clover.shape)
                
        Clover = preprocessing.scale(Clover,axis = 1)
        
        pca = PCA(n_components=3, whiten=False)
        COMPS = pca.fit_transform(Clover)
        X_se = pca.components_.T*np.sqrt(pca.explained_variance_)
        
        ###############################################################################
        ########################### PAINTING SHIT!! ###################################
        ###############################################################################
        
        ## 
        ## CLUSTER EIGENVALUES
        ##
        ### Clustering on decomposition
        from sklearn.cluster import KMeans
        kmeans = KMeans(n_clusters=9, random_state=0).fit(COMPS)
        labels1 = kmeans.labels_
        label_select = {y:[x for x in range(len(labels1)) if labels1[x] == y] for y in sorted(list(set(labels1)))}
        
        Cameo = []
        
        for cramp in sorted(label_select.keys()):
            Clamp = np.mean(Test[label_select[cramp],:],axis = 0)
            Fry = [Clamp[x] for x in range(len(Trend))]
            Cameo.append(Fry)
        
        Cameo = np.array(Cameo).T
        
        COMPS = pd.DataFrame(COMPS,columns=['PC1','PC2','PC3'])
        COMPS['order'] = [x for x in range(len(COMPS))]
        COMPS['labels'] = labels1
        X_se= pd.DataFrame(X_se,columns= ['PC1','PC2','PC3'])
        X_se['order'] = [x for x in range(len(X_se))]
        X_se['labels'] = Trend
        Cameo = pd.DataFrame(Cameo)
        Cameo['order'] = [x for x in range(len(Cameo))]
        
        return {
        'clusters': COMPS.to_json(orient='columns'),
        'loadings': X_se.to_json(orient = 'columns'),
        'Likes': Cameo.to_json(orient='columns')
        }



## read prepared ideogram:
file_name='example_chr1.pkl'
ideogram_bl = pd.read_pickle(file_name)

########

#### reading from test json (Much faster.. hopefully):

color_ref= ['red','yellow','blue','black','green','purple','orange','deepskyblue2','red3','darkolivegreen1','navy','chartreuse','darkorchid3','goldenrod2']
##
orderCore= pd.read_csv('Order_core_csv.txt')
Trend = [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]


### markdown

markdown_intro= '''
### Admixed Origin Dash Application

Explore the origin of hybrid accessions.

* The code for this application is available on [github](https://github.com/Joaos3092/Ideogram_vis)
\n

### Guide
Below you will find a description of targeted genetic variation at three specific
regions of chromosome 1 of *Oryza Sativa*. 
\n
Preceding this analysis, a whole genome crawl was performed to in order to assess the \n
most likely origin, in population terms, of each region of each accession in the data set. 
\n
The first graph, if `View` is set to `ALL`, is the output of that crawl along Chr 1 for 40 cBasmati accessions.
\n
the colors represent classifications into reference populations, allowing for 2 and 3-way uncertainty:

   
```
group_colors= { \n
    "blue": "Japonica" \n
    "yellow": "circumAus" \n
    "red": "Indica" \n
    "purple": "Indica-Japonica" \n
    "orange": "Indica-cAus" \n
    "green": "cAus-Japonica" \n
    "silver": "cAus-Indica-Japonica" \n
    "black": "outlier" \n
}
```


\n
This colorful plot is the first output of our exploration into the origin of these accessions. \n
However, we would also like to know if it is possible to identify subsets of the populations of origin closer \n
to the actual donors of this introgressed material.
\n

the question we try to answer here is:
**What structure underlies classification into a particular class?**\n

Using this application, you can explore the truth behind the output of the crawl and the classifications obtained.\n
After looking at the overall classifications, use the range slider below to select a region.\n
At the same time, chose what kind of classification you wish to focus on in the `Requested class` dropdown.\n

Profiles of the clusters each was connected with are then extracted. If `View` is set to `Requested`, then only the location of these\n
selected clusters is shown, in the colour of the chosen reference c. \n 
What follows is an analysis of the correlations among those clusters, and of what it tells us about genetic affiliations in the chosen regions.

**Attention:** Free heroku servers are limited, making the app lag. In case one of the figures doesn't update\n
immediately, just switch some of the parameters to get it to react.
'''


app = dash.Dash(__name__)

server = app.server

app.css.append_css({'external_url': 'https://codepen.io/chriddyp/pen/dZVMbK.css'})



app.layout = html.Div([
    
    html.Div([
    dcc.Markdown(children= markdown_intro)]),
    
    html.Hr(),
    
    html.Div([
    html.H3(
    id= "header1",className= "six columns",children= "Hello again"
    )],className= "row"),
    
    
    html.Div(
    id= "ideogram"
    ),
    
    html.Div([
    html.H6(children= 'View'),
    dcc.RadioItems(
    id= 'Analysis',
    className= 'four columns',
    value= 1,
    labelStyle={'display': 'inline-block'},
    options = [{'label':'All','value': 0},
               {'label':'Indica','value': 1},
                {'label':'Aus','value': 2},
                {'label':'Japonica','value': 3}]
    ),
    html.Button('THINK', id='button', n_clicks=0),
    html.Button('SHOW', id='button2', n_clicks=0)],className= 'row'),
    
    html.Div([
    dcc.RangeSlider(
    id= 'Region',
#    count= 1,
    min=0,
    max=44,
    step=1,
#    marks= {i:str(i) for i in np.arange(0,45,3)},
    value= [10,11]
    )
    ]),
    
    html.Hr(),
    
    
    html.Div([
    html.H4(
    id= "header1",className= "six columns",children= "Feature space"
    ),
    html.H4(
    id= "header2",className= "six columns",children= "selected accessions"
    )    
    ],className= "row"),
    
    
    html.Div([
    dcc.Graph(id='local_pca',animate = True,className= 'six columns'),
    html.Div(id='table-content', className='six columns')],className='row'),
    
    html.Div([
    html.Div(dcc.Markdown(children="""**Fig. 1** relative distances among accessions given cluster profiles selected and analysed. 
    In fact, loadings plot\n of PCA run on the former."""),className='six columns'),
    html.Div(dcc.Markdown(children= """**Table 1** Passport information on accessions shown in Fig. 1. If cluster cloud is selected
    below, then only accessions in red (updated plot) are shown."""),className= 'six columns')
    ],className= "row"),
    
    html.Hr(),
    
    html.Div([
    html.H5(
    id= "header1",className= "six columns",children= "opacity"
    ),
    html.H5(
    id= "header2",className= "six columns",children= "Likelihood threshold"
    )    
    ],className= "row"),
    
    html.Div([

    dcc.Slider(
    updatemode= "drag",
    className='six columns',
    id = 'opacity',
    min = .1,
    max = 1,
    value = .8,
    step = .1,
    #marks = [str(.1*x) for x in range(1,10)]
    ),
    dcc.Slider(
    updatemode= "drag",
    className= "six columns",
    id= "threshold",
    value= .1,
    min= .05,
    max= 1,
    marks= {i:str(i) for i in np.arange(.05,1,.05)},
    step= .05
    )
    ],className='row'),
    
    html.Hr(),
    
    html.Div([
    html.H5(children= 'Chose a color')
    ],className= "row"),
    
    html.Div(dcc.Markdown(children= """Clusters profiles in **Fig. 3** were grouped by colour and given a number each (hover over the points to see number).
    Proportion of different cluster types is plotted below to help you chose interesting clusters to analyse.
    Proximity among cluster profiles correlates with individual contribution patterns. When a group is chosen, the density of mean individual likelihoods across
    that group's profiles is plotted below (**Fig. 5**). Accessions with mean likelihoods above the *Lilelihood threshold* selected above will appear in red in **Fig. 1**."""),className= 'row'),
    
    
    html.Div([
    dcc.Dropdown(
    className='six columns',
    id = 'chose_color',
    value = 0,
    options = [{"label":x,"value": x} for x in range(10)]
    )
    ],className= "row"),
    
    html.Hr(),
    
    html.Div([
    dcc.Graph(id = "clusters",className="six columns"),
    dcc.Graph(id= "density_plot",className= "six columns")
    ]),
    
    html.Div([html.Div(dcc.Markdown(children= """**Fig.3** Relation among cluster profiles selected. In fact the distribution in feature space 
    of these profiles\n following principal component analysis."""),className= 'six columns'),
    html.Div(dcc.Markdown(children= """**Fig. 5** Density plot of average Likelihood by accession across cluster cloud selected."""),className= 'six columns')],
    className= 'row'),
    
    html.Div([
    dcc.Graph(id= "bars",className= "six columns"),
    
    html.Div(html.Div(dcc.Markdown(children= """**Fig. 4** proportion of cluster profiles by cloud (read *cluster*) in **clusters - observations**"""),
                      className= 'six columns'),className= 'row'),
    
    ]),
    
    
    html.Div(id= 'intermediate_labels',style= {'display': 'none'}),
    
])


@app.callback(
    Output('intermediate_labels','children'),
    [Input('Analysis','value'),
     Input('Region','value'),
    Input('button','n_clicks')]
)
def get_Clover(which,Region,button):
    if which != 0 and button >0:
        target= [which]
        Focus = ['CX59', 'CX65', 'CX67', 'CX104', 'CX143', 'CX149', 'IRIS_313-8268', 'IRIS_313-8326', 'IRIS_313-8385', 'IRIS_313-8656', 'IRIS_313-8712', 'IRIS_313-8747', 'IRIS_313-8813', 'IRIS_313-9083', 'IRIS_313-9172', 'IRIS_313-9601', 'IRIS_313-9629', 'IRIS_313-10670', 'IRIS_313-10851', 'IRIS_313-10868', 'IRIS_313-10926', 'IRIS_313-10933', 'IRIS_313-11021', 'IRIS_313-11022', 'IRIS_313-11023', 'IRIS_313-11026', 'IRIS_313-11218', 'IRIS_313-11258', 'IRIS_313-11268', 'IRIS_313-11289', 'IRIS_313-11293', 'IRIS_313-11451', 'IRIS_313-11564', 'IRIS_313-11567', 'IRIS_313-11625', 'IRIS_313-11627', 'IRIS_313-11630', 'IRIS_313-11632', 'IRIS_313-11743', 'IRIS_313-11825']
        chromosomes= [1]
        cons_threshold= .8
        
        Profiles,Names= db_profile_query('test_database.db',[x*1e6 for x in Region])
        
        focus_indexes= [x for x in range(len(Names)) if Names[x] in Focus]
        X_threshold= 5e-2
        
        Regard=  Clover(Profiles,focus_indexes,Trend,chromosomes,cons_threshold,X_threshold,Focus,target,Region)
        del Profiles
        return Regard


@app.callback(
    Output('ideogram','children'),
    [Input('Analysis','value'),
     Input('Region','value')]
)
def return_Ideogram(which,region):
    if which == 0:
        image_filename = 'Ideo_IRIS_313-11825.png'
        encoded_image = base64.b64encode(open(image_filename, 'rb').read())
        return [html.Img(id= 'spore',src='data:image/png;base64,{}'.format(encoded_image.decode()))]
    else:
        
        ideo = ideogram_bl.loc[(ideogram_bl.start >= region[0]*1e6) & (ideogram_bl.end < region[1]*1e6) & (ideogram_bl.gieStain == color_ref[which-1])]
                        
        layout = {'autosize': True, 'hovermode': 'closest', 'margin': {'l': 250, 'r': 199, 't': 120, 'b': 109, 'pad': 0}, 'xaxis1': {'anchor': 'y1', 'zeroline': False, 'ticks': 'inside', 'type': 'linear', 'range': [-2161775.8500000001, 45397292.850000001], 'showgrid': False, 'domain': [0.0, 1.0], 'side': 'bottom', 'tickfont': {'size': 10.0}, 'tick0': 0, 'dtick': 2000000, 'tickmode': False, 'mirror': 'ticks', 'showline': True}, 'yaxis1': {'anchor': 'x1', 'zeroline': False, 'ticks': 'inside', 'type': 'linear', 'range': [-1.0975000000000008, 23.047500000000014], 'showgrid': False, 'domain': [0.0, 1.0], 'side': 'left', 'tickfont': {'size': 10.0}, 'tick0': 21.700000000000014, 'dtick': -0.55000000000000071, 'tickmode': False, 'mirror': 'ticks', 'showline': True}, 'showlegend': False}
        #'width': 2000, 'height': 1000, 
        return [dcc.Graph(id= 'spore',figure = return_figure(ideo,layout))]


@app.callback(
    Output('bars','figure'),
    [Input('intermediate_labels','children'),
     Input('button2','n_clicks')])
def cluster_bars(Load,button):
    if Load and button > 0:
        clusters= pd.read_json(Load['clusters'])
        whom= sorted(list(set(clusters['labels'])))
        print('cluster bars')
        print(whom)
        nb= [round(len([x for x in clusters['labels'] if x == y]) / float(len(clusters)),3) for y in whom]
        nc= [str(x + 1) for x in whom]
        trace = [go.Bar(
        x= nc,
        y= nb,
        text= nb,
        marker=dict(
            color='rgb(158,202,225)',
            line=dict(
                color='rgb(8,48,107)',
                width=1.5),
        ),
        opacity= .6
        )]
        layout= go.Layout(
        title= 'cluster proportions'
        )
        fig= go.Figure(data=trace,layout=layout)
        return fig



@app.callback(
   Output('density_plot','figure'),
    [Input('chose_color','value'),
     Input('intermediate_labels','children'),
    Input('button2','n_clicks')])
def update_density(selected_group,Load,button):
    if Load and button > 0 and selected_group != 0:
        vectors = pd.read_json(Load['Likes'])
        vectors= vectors.sort_values("order")
        dense_plot=  ff.create_distplot([vectors.iloc[:,selected_group-1]], [str(selected_group)])
        dense_plot['layout'].update(title='<b>likelihood density</b>')
        return dense_plot



@app.callback(
    Output('local_pca','figure'),
    [Input('intermediate_labels','children'),
     Input('button2','n_clicks'),
    Input('threshold','value'),
    Input('opacity','value'),
    Input('chose_color','value')]
)
def update_loadings(Load,button,threshold,opac,selected_column):
    print('dangit ny')
    if Load and button > 0:
#        selected_column = 0
#        threshold = .2
#        opac = .8
        print('why')
        
        df= pd.read_json(Load['loadings'])
        df= df.sort_values("order")
        vectors= pd.read_json(Load['Likes'])
        vectors= vectors.sort_values("order")
        if selected_column == 0:
            scheme = Trend
            coords = {y:[x for x in range(len(scheme)) if scheme[x] == y] for y in list(set(scheme))}
            pop_refs= ["Indica","cAus","Japonica","GAP","cBasmati","Admix"]
            color_here= color_ref
        else:
            scheme = [int(vectors.iloc[x,selected_column-1]>=threshold) for x in range(len(df))]
            coords = {y:[x for x in range(len(scheme)) if scheme[x] == y] for y in list(set(scheme))}
            pop_refs= ["Below threshold","Above threshold"]
            color_here= ["grey","red"]
        
        figurine = loadings_graph(df,orderCore,scheme,coords,pop_refs,color_here,opac)
        return figurine


@app.callback(
    Output('table-content','children'),
    [Input('intermediate_labels','children'),
    Input('button2','n_clicks'),
    Input('threshold','value'),
    Input('chose_color','value')])
def update_table(Load,button,threshold,selected_group):
    if Load == None:
        show_table = [x for x in range(len(Trend))]
    else:
        selected_group = 1
        threshold = .2
        vectors= pd.read_json(Load['Likes'])
        vectors = vectors.sort_values("order")
        show_table = [x for x in range(len(vectors)) if vectors.iloc[x,selected_group-1] >= threshold]
    
    return [html.Div(
        id= 'table',
        children= generate_table(orderCore[["ID","NAME","COUNTRY","Initial_subpop"]].iloc[show_table,:]),
        style={
            'overflowX': 'scroll',
            'overflowY': 'scroll',
            'height': '450px',
            'display': 'block',
            'paddingLeft': '15px'
        })]


@app.callback(
    Output("clusters","figure"),
    [Input('intermediate_labels','children'),
     Input('button2','n_clicks')])
def update_secondFigure(Load,button):
    if Load and button > 0:
        cluster_pca= pd.read_json(Load['clusters'])
        return {'data': [go.Scatter3d(
            x = cluster_pca.iloc[[x for x in cluster_pca['labels'] if x == i]]['PC1'],
            y = cluster_pca.iloc[[x for x in cluster_pca['labels'] if x == i]]['PC2'],
            z = cluster_pca.iloc[[x for x in cluster_pca['labels'] if x == i]]['PC3'],
            type='scatter3d',
            mode= "markers",
            marker= {
    #            'color': [color_ref[x] for x in cluster_pca[0]],
    #            'color': cluster_pca[0],
                'line': {'width': 0},
                'size': 4,
                'symbol': 'circle',
                'opacity': .8
              },
              name = i + 1
            ) for i in cluster_pca['labels'].unique()],
            'layout': {
          "autosize": True, 
          "hovermode": "closest",
          "legend": {
            "x": 0.873529411765, 
            "y": 0.877829326396, 
            "borderwidth": 1, 
            "font": {"size": 13}
          },
          "scene": {
            "aspectmode": "auto", 
            "aspectratio": {
              "x": 1.02391505715, 
              "y": 0.737436541286, 
              "z": 1.3243763495
            }, 
            "camera": {
              "center": {
                "x": 0, 
                "y": 0, 
                "z": 0
              }, 
              "eye": {
                "x": 1.80578427889, 
                "y": 1.17729688569, 
                "z": 0.201532084509
              }, 
              "up": {
                "x": 0, 
                "y": 0, 
                "z": 1
              }
            }, 
            "xaxis": {
              "title": "PC1", 
              "type": "linear"
            }, 
            "yaxis": {
              "title": "PC2", 
              "type": "linear"
            }, 
            "zaxis": {
              "title": "PC3", 
              "type": "linear"
            }
          }, 
          "showlegend": False, 
          "title": "<b>clusters - observations</b>", 
          "xaxis": {"title": "V3"}, 
          "yaxis": {"title": "V2"}
        }
        }




# Run the Dash app
if __name__ == '__main__':
    app.server.run(debug=True,processes = 6)

