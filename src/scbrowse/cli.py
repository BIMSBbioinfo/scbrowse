"""
Module that contains the command line app.

Why does this file exist, and why not put this in __main__?

  You might be tempted to import things from __main__ later, but that will cause
  problems: the code will get executed twice:

  - When you run `python -mscbrowse` python will execute
    ``__main__.py`` as a script. That means there won't be any
    ``scbrowse.__main__`` in ``sys.modules``.
  - When you import __main__ it will get executed again (as a module) because
    there's no ``scbrowse.__main__`` in ``sys.modules``.

  Also see (1) from http://click.pocoo.org/5/setuptools/#setuptools-integration
"""
import io
import time
import base64
import os
from functools import wraps
from copy import copy
import logging
import dash  # pylint: disable=import-error
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import dash_core_components as dcc  # pylint: disable=import-error
import dash_html_components as html  # pylint: disable=import-error
import matplotlib.pyplot as plt
import plotly.graph_objs as go  # pylint: disable=import-error
import plotly.express as px
from flask_caching import Cache
from flask import Flask
import hashlib
import json
from multiprocessing import Pool
from scbrowse.utils import figure_encoder, draw_figure
import numpy as np
from pybedtools import BedTool

from anndata import read_h5ad
import scanpy as sc
import uuid
import coolbox
from coolbox.api import *
from coolbox.utilities import split_genome_range
from scbrowse.utils import SingleCellTrack

args = dict(matrix=os.environ['SCBROWSE_MATRIX'],
            genes=os.environ['SCBROWSE_GENES'],
            port=8051,
            logs=os.environ['SCBROWSE_LOGS'])

print(""" scbrowser startup .... """)
annotationcolors = sc.pl.palettes.vega_20_scanpy
selectioncolors = sc.pl.palettes.zeileis_28

###############
# data helper functions
###############

def log_layer(fn):
    @wraps(fn)
    def _log_wrapper(*args, **kwargs):
        ctx = dash.callback_context
        def _extr(props):
            return {k: props[k] for k in props \
                    if k not in ['points',
                             'scatter-plot.selectedData']}

        logging.debug(f'wrap:callb:'
                      f'outputs:{_extr(ctx.outputs_list)}:'
                      f'inputs:{_extr(ctx.inputs)}:'
                      f'triggered:{_extr(ctx.triggered[0])}:'
                     )
        try:
            return fn(*args, **kwargs)
        except PreventUpdate:
            raise
        except:
            logging.exception('callb:ERROR:')
            raise

    return _log_wrapper



logging.basicConfig(filename = args['logs'],
                    level=logging.DEBUG,
                    format='%(asctime)s;%(levelname)s;%(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

logging.debug('scbrowse - startup')
logging.debug(args)

##############
# load data
##############

genefile = args['genes']
genes = BedTool(args['genes'])

logging.debug(f'Number of genes: {len(genes)}')

ADATA = read_h5ad(args['matrix'])
ADATA.X.data[ADATA.X.data>0]=1


use_emb = ADATA.uns['embeddings'][0]
ADATA.var.loc[:, "nFrags"] = np.asarray(ADATA.X.sum(0)).flatten()
ADATA.obs.loc[:, "nFrags"] = np.asarray(ADATA.X.sum(1)).flatten()

ADATA.uns['nFrags'] = ADATA.var.nFrags.sum()


ADATA.obs.loc[:, "total"] = ADATA.obs.loc[:, "nFrags"]*1e5 / ADATA.uns['nFrags']
chroms = ADATA.obs.chrom.unique().tolist()

options = [dict(label=c, value=c) for c in chroms]
chromlens = {c: ADATA.obs.query(f'chrom == "{c}"').end.max() for c in chroms}

genelocus = [dict(label=g.name,
                  value=f'{g.chrom}:{max(1, g.start-10000)}-{min(chromlens[g.chrom], g.end+10000)}') for g in genes
             if g.chrom in chromlens]

for groups in ADATA.var.columns:
    if groups == 'nFrags':
        continue
    logging.debug(f'pre-compile {groups}')
    tracknames = sorted(ADATA.var[groups].unique().tolist())
    for i, track in enumerate(tracknames):
        sadata = ADATA[:,ADATA.var[groups]==track]
        ADATA.obs.loc[:, f'{groups}_{track}'] = \
                 np.asarray(sadata.X.sum(1)).flatten() *1e5 / sadata.var.nFrags.sum()

logging.debug(f'CountMatrix: {ADATA}')
print(repr(ADATA))
print(""" scbrowser running .... """)

##############
# instantiate app
##############

server = Flask("scbrowse")

app = dash.Dash("scbrowse", server=server)

app.title = "scbrowser"
app.config["suppress_callback_exceptions"] = True
cache = Cache(app.server, config={
    'CACHE_TYPE': 'redis',
    # Note that filesystem cache doesn't work on systems with ephemeral
    # filesystems like Heroku.
    'CACHE_TYPE': 'filesystem',
    'CACHE_DIR': 'cache-directory',

    # should be equal to maximum number of users on the app at a single time
    # higher numbers will store more data in the filesystem / redis cache
    'CACHE_THRESHOLD': 800
})

def get_cells(session_id, datahash, data=None):
    datahash = session_id + '_cells_' + datahash
    if data is not None:
        cache.set(datahash, data)
        logging.debug(f'get_cells({datahash}, <<<data>>>>)')
    logging.debug(f'get_cells({datahash})')
    return cache.get(datahash)

def get_selection(session_id, datahash, data=None):
    datahash = session_id + '_selection_' + datahash
    if data is not None:
        cache.set(datahash, data)
        logging.debug(f'get_selection({datahash}, <<<data>>>>)')
    logging.debug(f'get_selection({datahash})')
    return cache.get(datahash)



##############
# define layout
##############

def make_server():
    session_id = str(uuid.uuid4())
    logging.debug(f'new session_id = {session_id}')
    return html.Div(
    [
        html.Div(
            [
                html.Label(html.B("Annotation:")),
                dcc.Dropdown(
                    id="annotation-selector",
                    options=[
                        {"label": name, "value": name}
                        for name in ADATA.var.columns if name != 'nFrags'
                    ]
                    + [{"label": "None", "value": "None"}],
                    value="None",
                    style=dict(width="40%", display="inline-block"),
                ),
                html.Button('Clear all selections',
                    id="clear-selection",
                    n_clicks=0,
                    title='Clear all selections in the embedding',
                    style=dict(width="20%", display="inline-block"),
                ),
                html.Button('Undo last selection',
                    id="undo-last-selection",
                    n_clicks=0,
                    title='Clear last selections in the embedding',
                    style=dict(width="20%", display="inline-block"),
                ),
            ],
            style=dict(width="49%", display="inline-block"),
        ),
        html.Div(
            [
                html.Div([
                    html.Label(html.B("Gene:")),
                    dcc.Dropdown(id="gene-selector",
                             value=genelocus[0]['value'],
                             options=genelocus,
                             style=dict(width="70%", display="inline-block"),
                    )],
                ),
                html.Div([
                    html.Label(html.B("Locus:")),
                    dcc.Input(id="locus-selector",
                              type="text",
                              value="chr1:1-1000000",
                              style=dict(width="70%", display="inline-block"),),
                ]),
                html.Div([
                    html.Label(html.B("Zoom in:")),
                    html.Button("10x", id="zoom-in_10",
                              n_clicks=0,
                              style=dict(width="7%", display="inline-block"),),
                    html.Button("3x", id="zoom-in_3",
                              n_clicks=0,
                              style=dict(width="7%", display="inline-block"),),
                    html.Button("1.5x", id="zoom-in_15",
                              n_clicks=0,
                              style=dict(width="7%", display="inline-block"),),
                    html.Label(html.B("Zoom out:")),
                    html.Button("1.5x", id="zoom-out_15",
                              n_clicks=0,
                              style=dict(width="7%", display="inline-block"),),
                    html.Button("3x", id="zoom-out_3",
                              n_clicks=0,
                              style=dict(width="7%", display="inline-block"),),
                    html.Button("10x", id="zoom-out_10",
                              n_clicks=0,
                              style=dict(width="7%", display="inline-block"),),
                ]),
                html.Div([
                    html.Label(html.B("Move:")),
                    html.Button("<<<", id="move-left_3",
                              title="move left 95%",
                              n_clicks=0,
                              style=dict(width="10%", display="inline-block"),),
                    html.Button("<<", id="move-left_2",
                              title="move left 45%",
                              n_clicks=0,
                              style=dict(width="10%", display="inline-block"),),
                    html.Button("<", id="move-left_1",
                              title="move left 15%",
                              n_clicks=0,
                              style=dict(width="10%", display="inline-block"),),
                    html.Button(">", id="move-right_1",
                              title="move right 15%",
                              n_clicks=0,
                              style=dict(width="10%", display="inline-block"),),
                    html.Button(">>", id="move-right_2",
                              title="move right 45%",
                              n_clicks=0,
                              style=dict(width="10%", display="inline-block"),),
                    html.Button(">>>", id="move-right_3",
                              n_clicks=0,
                              title="move right 95%",
                              style=dict(width="10%", display="inline-block"),),
                ]),
            ],
            style=dict(width="45%", display="inline-block"),
        ),
        html.Br(),
        html.Div(
            [dcc.Graph(id="scatter-plot", style={'height': '500px'})],
            style=dict(width="49%", display="inline-block", verticalAlign='top'),
        id='divforscatter',),
        html.Div(
            [
                html.Img(id="genome-track", style=dict(verticalAlign='top')),
            ], id='divforgenome',
            style=dict(width="49%", display="inline-block"),
        ),
        html.Br(),
        html.Div([
                html.Div(id="stats-field"),
                dcc.Store(
                    id="session-id", data=session_id,
                ),
                dcc.Store(
                    id="dragmode-scatter", data="select",
                ),
                dcc.Store(
                    id="locus-store", data="chr1:1-1000000",
                ),
                dcc.Store(
                    id="selection-store", data=None,
                ),
        ],
        style=dict(width="49%", display="inline-block"),
        ),
    ],
    className="container",
)

app.layout = make_server
############
# app callbacks
############

def zoom(locus, factor):
    chrom, start, end = split_genome_range(locus)
    window = int((end - start)*factor/2)
    mid = (end+start)//2
    newl = f'{chrom}:{max(1, mid-window)}-{min(mid+window, chromlens[chrom])}'
    return newl


def move(locus, pshift):
    chrom, start, end = split_genome_range(locus)
    shift = int((end-start)*pshift)
    newl = f'{chrom}:{max(1, start+shift)}-{min(end+shift, chromlens[chrom])}'
    return newl

@app.callback(
    Output(component_id="locus-selector", component_property="value"),
    [
     Input(component_id="gene-selector", component_property="value"),
     Input(component_id="zoom-in_15", component_property="n_clicks"),
     Input(component_id="zoom-in_3", component_property="n_clicks"),
     Input(component_id="zoom-in_10", component_property="n_clicks"),
     Input(component_id="zoom-out_15", component_property="n_clicks"),
     Input(component_id="zoom-out_3", component_property="n_clicks"),
     Input(component_id="zoom-out_10", component_property="n_clicks"),
     Input(component_id="move-left_1", component_property="n_clicks"),
     Input(component_id="move-left_2", component_property="n_clicks"),
     Input(component_id="move-left_3", component_property="n_clicks"),
     Input(component_id="move-right_1", component_property="n_clicks"),
     Input(component_id="move-right_2", component_property="n_clicks"),
     Input(component_id="move-right_3", component_property="n_clicks"),
    ],
    [
     State(component_id="locus-store", component_property="data"),
    ],

)
@log_layer
def update_locus_selector_value(gcoord,zi15,zi5,zi10,zo15,zo5,zo10,ml1,ml2,ml3,mr1,mr2,mr3,storelocus):

    ctx = dash.callback_context
    if ctx.triggered is None:
        raise PreventUpdate

    if ctx.triggered[0]['prop_id'] == 'gene-selector.value':
        return gcoord

    if ctx.triggered[0]['prop_id'] == 'zoom-in_15.n_clicks':
        return zoom(storelocus, 1/1.5)
    if ctx.triggered[0]['prop_id'] == 'zoom-in_3.n_clicks':
        return zoom(storelocus, 1/3)
    if ctx.triggered[0]['prop_id'] == 'zoom-in_10.n_clicks':
        return zoom(storelocus, 1/10)

    if ctx.triggered[0]['prop_id'] == 'zoom-out_15.n_clicks':
        return zoom(storelocus, 1.5)
    if ctx.triggered[0]['prop_id'] == 'zoom-out_3.n_clicks':
        return zoom(storelocus, 3)
    if ctx.triggered[0]['prop_id'] == 'zoom-out_10.n_clicks':
        return zoom(storelocus, 10)

    if ctx.triggered[0]['prop_id'] == 'move-left_1.n_clicks':
        return move(storelocus, -.15)
    if ctx.triggered[0]['prop_id'] == 'move-left_2.n_clicks':
        return move(storelocus, -.45)
    if ctx.triggered[0]['prop_id'] == 'move-left_3.n_clicks':
        return move(storelocus, -.95)

    if ctx.triggered[0]['prop_id'] == 'move-right_1.n_clicks':
        return move(storelocus, .15)
    if ctx.triggered[0]['prop_id'] == 'move-right_2.n_clicks':
        return move(storelocus, .45)
    if ctx.triggered[0]['prop_id'] == 'move-right_3.n_clicks':
        return move(storelocus, .95)

    raise PreventUpdate



@app.callback(
    Output(component_id="scatter-plot", component_property="figure"),
    [
     Input(component_id="annotation-selector", component_property="value"),
     Input(component_id="selection-store", component_property="data"),
    ],
    [
     State(component_id="dragmode-scatter", component_property="data"),
     State(component_id="session-id", component_property="data"),
    ]
)
@log_layer
def embedding_callback(annot, selection_store, dragmode, session_id):
    co = {}
    if annot in ADATA.var.columns:
        tracknames = sorted(ADATA.var[annot].unique().tolist())
    else:
        tracknames = ['None']

    colors = {trackname: annotationcolors[i%len(annotationcolors)] for i, \
              trackname in enumerate(tracknames)}
    colors['None'] = 'black'

    if selection_store is not None and len(selection_store) > 0:
       selection_hash = selection_store[-1]
       selection = get_selection(session_id, selection_hash)
       selnames = [selkey for selkey in selection]
       colors.update({trackname: selectioncolors[i%len(selectioncolors)] for i, \
                 trackname in enumerate(selnames)})
       tracknames += selnames

    df = ADATA.var.copy()
    df.loc[:, 'dim1'] = ADATA.varm[use_emb][:,0]
    df.loc[:, 'dim2'] = ADATA.varm[use_emb][:,1]

    fig = px.scatter(
        df,
        x="dim1",
        y="dim2",
        hover_name=df.index,
        opacity=0.3,
        color=annot if annot != "None" else None,
        color_discrete_sequence=['black'] if annot == "None" else None,
        color_discrete_map=colors,
        custom_data=[df.index],
        template="plotly_white",
    )

    if selection_store is not None and len(selection_store) > 0:
       selection_hash = selection_store[-1]
       selection = get_selection(session_id, selection_hash)

       for selkey in selection:
           fig.add_trace(go.Scatter(x=selection[selkey]['x'],
                                    y=selection[selkey]['y'],
                                    text=selkey,
                                    name=selkey,
                                    fillcolor=colors[selkey],
                                    opacity=0.2,
                                    mode='lines',
                                    hoverinfo='skip',
                                    line_width=0,
                                    fill='toself'))
    fig.update_traces(marker=dict(size=3))
    fig.update_layout(clickmode="event+select",
                      template='plotly_white',
                      dragmode=dragmode)

    return fig

pool = Pool(10)

@app.callback(
    Output(component_id="genome-track", component_property="src"),
    [
        Input(component_id="locus-selector", component_property="value"),
        Input(component_id="annotation-selector", component_property="value"),
        Input(component_id="selection-store", component_property="data"),
    ],
    [
        State(component_id="session-id", component_property="data"),
    ],
)
@log_layer
def genome_tracks_callback(locus,
                   #highlight, overlay,
                   annotation, selectionstore, session_id):
    if locus is None:
        raise PreventUpdate

    if selectionstore is not None:
        selcells = get_cells(session_id, selectionstore[-1])
    else:
        selcells = None

    chrom, start, end = split_genome_range(locus)

    sada = ADATA[(ADATA.obs.chrom==chrom) & (ADATA.obs.start>=start) & (ADATA.obs.end<=end),:].copy()
    #return draw_figure((ADATA, locus, annotation, selcells, genefile))

    ts = time.time()
    print('using pool')
    ret = pool.map(draw_figure, ((sada, locus, annotation, selcells, genefile),))
    print(f'pool {time.time()-ts}')
    return ret[0]
    #print(ret)


@app.callback(
    Output(component_id="locus-store", component_property="data"),
    [
     Input(component_id="locus-selector", component_property="value"),
    ],
)
@log_layer
def locus_store_callback(coord):

    ctx = dash.callback_context
    if ctx.triggered is None:
        raise PreventUpdate

    if  ctx.triggered[0]['prop_id'] == 'locus-selector.value':

        return coord

    raise PreventUpdate


@app.callback(
    Output(component_id="selection-store", component_property="data"),
    [
        Input(component_id="scatter-plot", component_property="selectedData"),
        Input(component_id="clear-selection", component_property="n_clicks"),
        Input(component_id="undo-last-selection", component_property="n_clicks"),
    ],
    [
        State("selection-store", "data"),
        State(component_id="session-id", component_property="data"),
    ],
)
@log_layer
def selection_store(selected, clicked, undolast, prev_hash, session_id):

    ctx = dash.callback_context

    if ctx.triggered is None:
        raise PreventUpdate

    if ctx.triggered[0]['prop_id'] == 'clear-selection.n_clicks':
        new_hash = [""]
        get_cells(session_id, new_hash[-1])
        get_selection(session_id, new_hash[-1], dict())
        return new_hash

    if ctx.triggered[0]['prop_id'] == 'undo-last-selection.n_clicks':
        if prev_hash is not None and len(prev_hash) > 1:
            new_hash = prev_hash[:-1]
        else:
            # don't pop if already empty
            new_hash = [""]
            get_cells(session_id, new_hash[-1])
            get_selection(session_id, new_hash[-1], dict())
        return new_hash

    if selected is None:
        raise PreventUpdate

    if prev_hash is not None:
        prev_cells = copy(get_cells(session_id, prev_hash[-1])) or {}
        prev_selection = copy(get_selection(session_id, prev_hash[-1])) or {}
    else:
        prev_hash = []
        prev_cells = dict()
        prev_selection = dict()

    selname = f'sel_{len(prev_cells)}'
    # got some new selected points
    sel = [point["customdata"][0] for point in selected["points"]]
    cell_ids = {selname:
                ADATA.var[ADATA.var.index.isin(sel)].index.tolist()}

    selectionpoints = dict()
    if 'range' in selected:
        name = 'range'
        xs = [point for point in selected[name]['x']]
        ys = [point for point in selected[name]['y']]
        selectionpoints[selname] = dict(x=[xs[0], xs[0], xs[1], xs[1], xs[0]],
                                        y=[ys[0], ys[1], ys[1], ys[0], ys[0]])
    if 'lassoPoints' in selected:
        name = 'lassoPoints'
        selectionpoints[selname] = \
           {k: [point for point in selected[name][k]] for k in selected[name]}
        selectionpoints[selname]['x'].append(selectionpoints[selname]['x'][0])
        selectionpoints[selname]['y'].append(selectionpoints[selname]['y'][0])

    new_cells = copy(prev_cells)
    new_cells.update(cell_ids)
    new_selection = copy(prev_selection)
    new_selection.update(selectionpoints)

    data_md5 = hashlib.md5(str(json.dumps(
        new_cells, sort_keys=True)).encode('utf-8')).hexdigest()

    new_hash = prev_hash + [data_md5]

    get_cells(session_id, new_hash[-1], new_cells)

    get_selection(session_id, new_hash[-1], new_selection)

    return new_hash


@app.callback(
    Output(component_id="dragmode-scatter", component_property="data"),
    [Input(component_id="scatter-plot", component_property="relayoutData"),
     ],
)
@log_layer
def update_dragmode_scatter(relayout):
    if relayout is None:
        raise PreventUpdate
    if "dragmode" not in relayout:
        raise PreventUpdate
    return relayout["dragmode"]


def main():
    ############
    # run server
    ############
    app.run_server(debug=True, port=args['port'])

if __name__ == '__main__':
    main()
