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
from collections import Counter
from copy import copy
import logging
import dash  # pylint: disable=import-error
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import dash_core_components as dcc  # pylint: disable=import-error
import dash_html_components as html  # pylint: disable=import-error
import plotly.graph_objs as go  # pylint: disable=import-error
import plotly.express as px

from plotly.subplots import make_subplots
from flask_caching import Cache
from flask import Flask
import hashlib
import json
import numpy as np
from pybedtools import BedTool, Interval

from anndata import read_h5ad
import scanpy as sc
import uuid

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

def split_genome_range(gr):
    chrom, rest = gr.split(':')
    start, end = rest.split('-')
    return chrom, int(start), int(end)

logging.basicConfig(filename = args['logs'],
                    level=logging.DEBUG,
                    format='%(asctime)s;%(levelname)s;%(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

logging.debug('scbrowse - startup')
logging.debug(args)

def load_dataset(h5file):
    adata = read_h5ad(h5file)
    adata.X.data[adata.X.data>0]=1

    adata.var.loc[:, "nFrags"] = np.asarray(adata.X.sum(0)).flatten()
    adata.obs.loc[:, "nFrags"] = np.asarray(adata.X.sum(1)).flatten()

    adata.uns['nFrags'] = adata.var.nFrags.sum()


    adata.obs.loc[:, "total"] = adata.obs.loc[:, "nFrags"]*1e5 / adata.uns['nFrags']

    for groups in adata.var.columns:
        if groups == 'nFrags':
            continue
        logging.debug(f'pre-compile {groups}')
        tracknames = sorted(adata.var[groups].unique().tolist())
        for i, track in enumerate(tracknames):
            sadata = adata[:,adata.var[groups]==track]
            adata.obs.loc[:, f'{groups}_{track}'] = \
                     np.asarray(sadata.X.sum(1)).flatten() *1e5 / sadata.var.nFrags.sum()

    logging.debug(f'CountMatrix: {adata}')
    #print(repr(adata))
    return adata


class TrackManager:
    def __init__(self, adata, locus, annotation, selected_cells, genes):
        self.obs = adata.obs.copy()
        self.tracknames = ["total"]
        self.colnames = ["total"]
        self.colors = {"total":"black"}

        if annotation != "None":
            #self.annotation = annotation
            #names = sorted(adata.var[annotation].unique().tolist())
            names = TRACKNAMES[annotation]
            colors = annotationcolors[:len(names)]
            #self.tracknames = sorted(adata.var[annotation].unique())
            for i, name in enumerate(names):
                self.tracknames.append(name)
                self.colors[name] = colors[i]
                self.colnames.append(f'{annotation}_{name}')

        if selected_cells is not None:
            colors = selectioncolors[:len(selected_cells)]
            for i, name in enumerate(selected_cells):
                sadata = adata[:,adata.var.index.isin(selected_cells[name])]
                da = np.asarray(sadata.X.sum(1)).flatten() *1e5 / sadata.var.nFrags.sum()
                self.obs.loc[:,name] = da
                self.tracknames.append(name)
                self.colors[name] = colors[i]
                self.colnames.append(f'{name}')


        self.trackheight = 3
        self.chrom, self.start, self.end = split_genome_range(locus)
        self.locus = [self.chrom, self.start,self.end]
        self.genes = genes

        self.init_trace()

    def allocate(self):

        ntracks = len(self.tracknames)
        specs=[]
        for n in range(ntracks):
            specs += [[{"rowspan": self.trackheight}]]
            specs += [[None]]*(self.trackheight - 1)
        specs += [[{"rowspan": 10}]]
        specs += [[None]]*(10 - 1)

        fig = make_subplots(
            rows=ntracks * self.trackheight + 10,
            cols=1,
            specs=specs,
            shared_xaxes=True,
            vertical_spacing=0.00,
        )
        fig.update_xaxes(range=[self.start, self.end])
        return fig

    def draw_track(self):
        plobjs = []
        sregs_ = self.obs
        plobjs.append(go.Scatter(x=sregs_.start,
                                 y=sregs_[self.colnames[self.itrack]],
                                 line=dict(width=0.1, color=self.colors[self.tracknames[self.itrack]]),
                                 fill='tozeroy',
                                 mode="lines", name=self.tracknames[self.itrack],))
        self.ymax = max(self.ymax, sregs_[self.colnames[self.itrack]].max())

        return plobjs

    def init_trace(self):
        self.itrace = 1
        self.itrack = 0
        self.ymax = 0.0

    def next_trace(self):
        self.itrace += self.trackheight
        self.itrack += 1

    def extend_trace(self, fig, plobjs):
        for plobj in plobjs:
            fig.add_trace(plobj, row=self.itrace, col=1)
        return fig

    def draw_tracks(self, fig):
        for i, trackname in enumerate(self.tracknames):
            self.extend_trace(fig, self.draw_track())
            self.next_trace()

        return fig

    def draw(self):
        fig = self.allocate()
        self.init_trace()
        fig = self.draw_tracks(fig)

        title = f"Window: {self.locus[0]}:{self.locus[1]}-{self.locus[2]}"
        fig.layout.update(
            dict(title=title,
                 clickmode='none',
                 template="plotly_white")
        )
        for prop in fig.layout:
            if 'yaxis' in prop:
                fig.layout[prop]["range"] = [0,self.ymax]

        fig = self.draw_annotation(fig)
        fig.layout[f'yaxis{self.itrack+1}']["range"] = [0,4.]

        for prop in fig.layout:
            if 'xaxis' in prop:
                fig.layout[prop]["showgrid"] = False
                fig.layout[prop]["zeroline"] = False
                fig.layout[prop]["fixedrange"] = True
                fig.layout[prop]["range"] = [self.start, self.end]
                fig.layout[prop]["showticklabels"] = False
            if 'yaxis' in prop:
                fig.layout[prop]["showgrid"] = False
                fig.layout[prop]["zeroline"] = False
                fig.layout[prop]["fixedrange"] = True
        fig.layout[f"xaxis"]["side"] = "top"
        fig.layout[f"xaxis"]["showticklabels"] = True
        return fig

    def draw_annotation(self, fig):
        ntracks = len(self.tracknames)
        chrom, start, end = self.locus
        plobjs = _draw_gene_annotation(fig, self.genes, chrom, start, end)
        for plobj in plobjs or []:
            fig.add_trace(
                plobj,
                row=ntracks*self.trackheight + 1,
                col=1,
            )
        if plobjs is not None:
            fig.layout[f"yaxis{ntracks+1}"]["showticklabels"] = False
            fig.layout[f"yaxis{ntracks+1}"]["showgrid"] = False
            fig.layout[f"yaxis{ntracks+1}"]["zeroline"] = False
            fig.layout[f"xaxis{ntracks+1}"]["showgrid"] = False
            fig.layout[f"xaxis{ntracks+1}"]["zeroline"] = False
        return fig


def draw_box(yoffset, start, end, height=0.1):
    xs = [
        start,
        start,
        end,
        end,
        start,
        None,
    ]
    ys = [yoffset, yoffset+height, yoffset+height, yoffset, yoffset, None]
    return xs, ys

def draw_right_arrow_box(yoffset, start, end, height=0.1):
    xs = [
        start,
        start,
        max(start, end-5000),
        end,
        max(start,end-5000),
        start,
        None,
    ]
    ys = [yoffset, yoffset+height, yoffset+height, yoffset + height*.5, yoffset, yoffset, None]
    return xs, ys

def draw_left_arrow_box(yoffset, start, end, height=0.1):
    xs = [
        start,
        min(start + 5000, end),
        end,
        end,
        min(start + 5000, end),
        start,
        None,
    ]
    ys = [yoffset +height*.5, yoffset, yoffset, yoffset + height, yoffset +height, yoffset + height*.5, None]
    return xs, ys


def draw_gene(yoffset, gene):
    return draw_simple_gene(yoffset, gene)

def draw_simple_gene(yoffset, gene):
    if gene.strand == "-":
        return draw_left_arrow_box(yoffset, gene.start, gene.end)
    else:
        return draw_right_arrow_box(yoffset, gene.start, gene.end)

def _draw_gene_annotation(fig, genes, chrom, start, end):
    wbed = BedTool([Interval(chrom, start, end)])
    regions = genes.intersect(wbed, wa=True, u=True)
    xs = []
    ys = []

    textxpos = []
    textypos = []
    names = []
    offset = 3.5-.07
    prevends = 0

    rangeannot = []
    for i, region in enumerate(regions):
        names.append(region.name)

        #if region.start >= (prevends +10000):
        #    offset = 3.5-.07

        x, y = draw_gene(offset, region)
        xs +=x
        ys +=y
        textxpos.append(region.end+500)
        textypos.append(offset)
        rangeannot.append(f"{region.chrom}:{region.start}-{region.end};{region.strand}")

        prevends = max(prevends, region.end)
        offset -= 0.2
        if offset <= 0.0:
            offset=3.5-0.07


    if len(textxpos) > 0:
        plobjs = [
            go.Scatter(
                x=xs,
                y=ys,
                mode="lines",
                fill="toself",
                name="Genes",
                marker=dict(color="goldenrod"),
            ),
            go.Scatter(
                x=textxpos,
                y=textypos,
                text=names,
                mode="text",
                #opacity=0.0,
                name="Genes",
                customdata=rangeannot,
                hovertemplate="%{text}<br>%{customdata}",
                showlegend=False,
            ),
            #go.Scatter(
            #    x=xs,
            #    y=ys,
            #    mode="lines",
            #    #fill="toself",
            #    #name="Genes",
            #    #marker=dict(color="black"),
            #    line=dict(color='black', width=.1),
            #    showlegend=False,
            #),
        ]
        return plobjs

##############
# load data
##############

genefile = args['genes']
genes = BedTool(args['genes'])

logging.debug(f'Number of genes: {len(genes)}')

ADATA = load_dataset(args['matrix'])
use_emb = ADATA.uns['embeddings'][0]
chroms = ADATA.obs.chrom.unique().tolist()

chromlens = {c: ADATA.obs.query(f'chrom == "{c}"').end.max() for c in chroms}
ANNOTATION = [name for name in ADATA.var.columns if name != 'nFrags']

TRACKNAMES = {}

for annot in ANNOTATION:
    TRACKNAMES[annot] = sorted(ADATA.var[annot].unique().tolist())

TRACKNAMES['None'] = ['None']
ADATA_subs = {chrom: ADATA[ADATA.obs.chrom==chrom,:] for chrom in chroms}

genelocus = [dict(label=g.name,
                  value=f'{g.chrom}:{max(1, g.start-10000)}-{min(chromlens[g.chrom], g.end+10000)}') for g in genes
             if g.chrom in chromlens]

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
                        for name in ANNOTATION
                    ]
                    + [{"label": "None", "value": "None"}],
                    value="cell_label",
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
            [
             dcc.Graph(id="scatter-plot", style={'height': '500px'}),
             html.P("Tip: Select cells in the embedding manually using the box or lasso selection tool."),
             html.P(["Source code: ",
                     html.A("scbrowse",
                     href="https://github.com/BIMSBbioinfo/scbrowse"),
                    ]
                   ),
            ],
            style=dict(width="49%",
                       display="inline-block",
                       verticalAlign='top'),
            id='divforscatter',),
        html.Div(
            [
                dcc.Graph(id="ply-genome-track", style={'height': '500px'},
                          config={'displayModeBar': False})
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
    tracknames = TRACKNAMES[annot]

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

@app.callback(
    Output(component_id="ply-genome-track", component_property="figure"),
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
def ply_genome_tracks_callback(locus,
                               annotation, selectionstore, session_id):
    if locus is None:
        raise PreventUpdate

    if selectionstore is not None:
        selected_cells = get_cells(session_id, selectionstore[-1])
    else:
        selected_cells = None

    st=time.time()
    chrom, start, end = split_genome_range(locus)

    sada = ADATA_subs[chrom]
    sada = sada[(sada.obs.start>=start) & (sada.obs.end<=end),:].copy()

    #obs = sada.obs
    tracks = TrackManager(sada, locus, annotation, selected_cells, genes)
    fig = tracks.draw()
    return fig


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
    app.run_server(debug=False, port=args['port'])

if __name__ == '__main__':
    main()
