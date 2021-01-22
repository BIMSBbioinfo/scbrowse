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
import os
import argparse
from functools import wraps
from copy import copy
import logging
import argparse
import dash  # pylint: disable=import-error
import dash_table
from dash.exceptions import PreventUpdate
import dash_core_components as dcc  # pylint: disable=import-error
import dash_html_components as html  # pylint: disable=import-error
import plotly.graph_objs as go  # pylint: disable=import-error
import plotly.express as px
from plotly.subplots import make_subplots
from scipy.stats import fisher_exact
from flask_caching import Cache
from flask import Flask
import hashlib
import json

from scipy.sparse import diags
import numpy as np
from dash.dependencies import Input, Output, State
from pybedtools import BedTool, Interval

import pandas as pd
from scregseg.countmatrix import CountMatrix
import uuid

args = dict(embedding=os.environ['SCBROWSE_EMBEDDING'],
            matrix=os.environ['SCBROWSE_MATRIX'],
            regions=os.environ['SCBROWSE_REGIONS'],
            genes=os.environ['SCBROWSE_GENES'],
            port=8051,
            logs=os.environ['SCBROWSE_LOGS'])

annotationcolors = px.colors.qualitative.Light24
selectioncolors = px.colors.qualitative.Dark24

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
                             'scatter-plot.selectedData',
                             'summary-plot.selectedData']}

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


def cell_coverage(cmat, regions, cells):
    m = cmat.cmat
    if regions is not None:
        m = m[regions, :]
    if cells is not None:
        #m = m[:, cmat.cannot[cmat.cannot.cell.isin(cells)].index]
        m = m[:, cells]
        return np.asarray(m.sum(1)).flatten()
    return np.asarray(m.sum(1)).flatten()


def _highlight(row, rmin, rmax, highlight):
    if (
        row.start >= rmin + (rmax - rmin) / 100 * highlight[0]
        and row.end <= rmin + (rmax - rmin) / 100 * highlight[1]
    ):
        return "inside"
    return "outside"

##############
# drawing helper functions
##############

def get_locus(regs_):
    if regs_.shape[0] == 0:
        return None
    xmin = regs_.start.min()
    xmax = regs_.end.max()
    chrom = regs_.chrom.unique()[0]
    return [chrom, xmin, xmax]

def get_highlight(regs_):
    sregs_ = regs_[regs_.highlight == "inside"]
    if sregs_.shape[0] == 0:
        return None
    xmin = sregs_.start.min()
    xmax = sregs_.end.max()
    return [xmin, xmax]


def data_bars_diverging(data, column, color_above='#3D9970', color_below='#FF4136'):
    data = [float(d) for d in data[column]]
    n_bins = 100
    bounds = [i * (1.0 / n_bins) for i in range(n_bins + 1)]
    col_max = 3.
    col_min = -3.
    ranges = [
        ((col_max - col_min) * i) + col_min
        for i in bounds
    ]
    midpoint = (col_max + col_min) / 2.

    styles = []
    for i in range(1, len(bounds)):
        min_bound = ranges[i - 1]
        max_bound = ranges[i]
        min_bound_percentage = bounds[i - 1] * 100
        max_bound_percentage = bounds[i] * 100

        style = {
            'if': {
                'filter_query': (
                    '{{{column}}} >= {min_bound}' +
                    (' && {{{column}}} < {max_bound}' if (i < len(bounds) - 1) else '')
                ).format(column=column, min_bound=min_bound, max_bound=max_bound),
                'column_id': column
            },
            'paddingBottom': 2,
            'paddingTop': 2
        }
        if max_bound > midpoint:
            background = (
                """
                    linear-gradient(90deg,
                    white 0%,
                    white 50%,
                    {color_above} 50%,
                    {color_above} {max_bound_percentage}%,
                    white {max_bound_percentage}%,
                    white 100%)
                """.format(
                    max_bound_percentage=max_bound_percentage,
                    color_above=color_above
                )
            )
        else:
            background = (
                """
                    linear-gradient(90deg,
                    white 0%,
                    white {min_bound_percentage}%,
                    {color_below} {min_bound_percentage}%,
                    {color_below} 50%,
                    white 50%,
                    white 100%)
                """.format(
                    min_bound_percentage=min_bound_percentage,
                    color_below=color_below
                )
            )
        style['background'] = background
        styles.append(style)

    return styles

class TableManager:
    def __init__(self, regs_):
        self.regs_ = regs_
        self.header = ["Name", "Highlight", "Outside", "LogOdds", "Pvalue"]
        self.data = {h: [] for h in self.header}
        self.nrows = 0

    def draw(self):
        if self.nrows <=0:
            return ""
        return dash_table.DataTable(
          id='datatable',
          columns=self.draw_header(),
          data=self.draw_rows(),
          style_cell={
             'whiteSpace': 'normal',
             'height': 'auto',
         },
          style_data_conditional=
          ( data_bars_diverging(self.data, 'LogOdds',
                                color_above=px.colors.diverging.Picnic[-1],
                                color_below=px.colors.diverging.Picnic[0])
          ),
          page_action='none',
          style_table={'height': '300px', 'overflowY': 'auto',
                       }
        )

    def draw_header(self):
        return [{'name': h, 'id': h, 'type': 'text' \
                 if h == 'Name' else 'numeric'} for h in self.header]

    def add_row(self, name, n11, c1, r1, n, ncells):
        n21 = c1 - n11
        n12 = r1 - n11
        n22 = n - n11 - n12 - n21
        odds, pval = fisher_exact([[n11, n12],
                                   [n21, n22]],
                                  alternative='greater')
        self.data[self.header[0]].append(f'{name} ({ncells})')
        self.data[self.header[1]].append(n11)
        self.data[self.header[2]].append(n21)
        self.data[self.header[3]].append('{:.3f}'.format(np.log2(odds)))
        self.data[self.header[4]].append('{:.3f}'.format(pval))
        self.nrows += 1

    def draw_rows(self):
        return [{h: self.data[h][irow] for h in self.header} \
                for irow in range(self.nrows)]



class TrackManager:
    def __init__(self, regs_, tracknames,
                 colors,
                 genes, controlprops):
        self.regs_ = regs_
        self.tracknames = tracknames
        self.trackheight = 3
        self.locus = get_locus(regs_)
        self.highlight = get_highlight(regs_)
        self.controlprops = controlprops
        self.genes = genes
        self.colors = colors
        assert len(tracknames) == len(colors)
        self.init_trace()

    def allocate(self):
        if not self.controlprops['overlay']:
            ntracks = 1
        else:
            ntracks = len(self.tracknames)
        specs=[]
        for n in range(ntracks):
            specs += [[{"rowspan": self.trackheight}]]
            specs += [[None]]*(self.trackheight - 1)
        specs += [[{"rowspan": 1}]]

        fig = make_subplots(
            rows=ntracks * self.trackheight + 1,
            cols=1,
            specs=specs,
            shared_xaxes=True,
            vertical_spacing=0.00,
        )
        return fig

    def draw_summary_track(self, trackname):
        plobjs = []
        plottype = self.controlprops['plottype']
        sregs_ = self.regs_
        if plottype == "bar":
            plobjs.append(go.Bar(x=sregs_.start,
                                 y=sregs_[trackname],
                                 marker=dict(color=self.colors[trackname]),
                                 name=trackname))
        else:
            plobjs.append(go.Scatter(x=sregs_.start,
                                     y=sregs_[trackname],
                                     marker=dict(color=self.colors[trackname]),
                                     mode="lines", name=trackname,))

        rangeannot = self.regs_.apply(lambda row: \
            f'{row.chrom}:{row.start}-{row.end}', axis=1).values.tolist()
        plobjs.append(go.Scatter(
                x=sregs_.start,
                y=sregs_[trackname],
                mode="markers",
                opacity=0,
                customdata=rangeannot,
                hovertemplate="y=%{y}<br>%{customdata}<br>" + trackname,
                showlegend=False,
            ))
        return plobjs

    def init_trace(self):
        self.itrace = 1

    def next_trace(self):
        if self.controlprops['overlay']:
            self.itrace += self.trackheight

    def extend_trace(self, fig, plobjs):
        for plobj in plobjs:
            fig.add_trace(plobj, row=self.itrace, col=1)
        return fig

    def draw_tracks(self, fig):
        self.init_trace()
        for trackname in self.tracknames:
            self.extend_trace(fig, self.draw_summary_track(trackname))
            self.next_trace()

        title = f"Window: {self.locus[0]}:{self.locus[1]}-{self.locus[2]}"
        if self.highlight is not None:
            title += f"<br>Highlight: {self.locus[0]}:{self.highlight[0]}-{self.highlight[1]}"

        fig.layout.update(
            dict(title=title, dragmode=self.controlprops['dragmode'],
                 clickmode="event+select",
                 template="plotly_white")
        )
        return fig

    def draw_highlight(self, fig):
        if self.highlight is None:
            return fig

        xmin, xmax = self.highlight
        return fig.add_shape(
            type="rect",
            xref='x',
            yref='paper',
            x0=xmin,
            y0=0,
            x1=xmax,
            y1=1,
            fillcolor="LightSalmon",
            opacity=0.3,
            layer="below",
            line_width=0,
        )

    def draw_annotation(self, fig):
        ntracks = len(self.tracknames) if self.controlprops['overlay'] else 1
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
        else:
            fig.layout["xaxis"]["showticklabels"] = True
        return fig

    def draw(self):
        fig = self.allocate()
        fig = self.draw_tracks(fig)
        fig = self.draw_highlight(fig)
        fig = self.draw_annotation(fig)

        return fig



def _draw_gene_annotation(fig, genes, chrom, start, end):
    wbed = BedTool([Interval(chrom, start, end)])
    regions = genes.intersect(wbed)
    xs = []
    ys = []

    midpoints = []
    names = []
    offset = 0
    lastiv = []

    rangeannot = []
    for i, region in enumerate(regions):
        names.append(region.name)

        # draw arrow to indicate direction
        if region.strand != "-":
            xs += [
                region.start,
                region.start,
                region.end,
                region.end + 1000,
                region.end,
                region.start,
                None,
            ]
            ys += [0, 1, 1, 0.5, 0, 0, None]
            midpoints.append(region.start)
        else:
            xs += [
                region.start,
                region.start - 1000,
                region.start,
                region.end,
                region.end,
                region.start,
                None,
            ]
            ys += [0, 0.5, 1, 1, 0, 0, None]
            midpoints.append(region.end)
        rangeannot.append(f"{region.chrom}:{region.start}-{region.end};{region.strand}")

    if len(midpoints) > 0:
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
                x=midpoints,
                y=[0.5] * len(midpoints),
                text=names,
                mode="text",
                opacity=0.0,
                name="Genes",
                customdata=rangeannot,
                hovertemplate="%{text}<br>%{customdata}",
                showlegend=False,
            ),
        ]
        return plobjs




logging.basicConfig(filename = args['logs'],
                    level=logging.DEBUG,
                    format='%(asctime)s;%(levelname)s;%(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

logging.debug('scbrowse - startup')
logging.debug(args)

##############
# load data
##############

genes = BedTool(args['genes'])

logging.debug(f'Number of genes: {len(genes)}')

emb = pd.read_csv(args['embedding'], sep="\t")
if "barcode" not in emb.columns:
    emb["barcode"] = emb.index

logging.debug(f'Embedding dims: {emb.shape}')

cmat = CountMatrix.from_mtx(args['matrix'], args['regions'])

readsincell = np.asarray(cmat.cmat.sum(0)).flatten()
readsinregion = np.asarray(cmat.cmat.sum(1)).flatten()

totcellcount = np.asarray(cmat.cmat.sum(0)).flatten()
totcellcount = totcellcount.astype("float") / totcellcount.sum()

logging.debug(f'CountMatrix: {cmat}')

regs = cmat.regions.copy()
regs["total_coverage"] = cell_coverage(cmat, None, None)

chroms = cmat.regions.chrom.unique().tolist()

options = [dict(label=c, value=c) for c in chroms]
chromlens = {c: regs.query(f'chrom == "{c}"').end.max() for c in chroms}

celldepth = np.asarray(cmat.cmat.sum(0)).flatten()

genelocus = [dict(label=g.name,
                  value=f'{g.chrom}:{g.start}-{g.end}') for g in genes]

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

def get_region_selection(chrom, interval, highlight):
    @cache.memoize()
    def _get_region_selection(chrom, interval, highlight):

        regs_ = regs.copy()
        start, end = interval

        regs_ = regs_.query(f'chrom == "{chrom}" & start >= {int(start) -1}')

        end = max(regs_.iloc[0, :].end, end)
        regs_ = regs_.query(f'end <= {end}')
        regs_.insert(regs_.shape[1], "annot", "other", True)

        rmin, rmax = regs_.start.min(), regs_.end.max()

        regs_.insert(
            regs_.shape[1],
            "highlight",
            regs_.apply(_highlight, axis=1, args=(rmin, rmax, highlight)),
            True,
        )
        return regs_.to_json()

    return pd.read_json(_get_region_selection(chrom, interval, highlight))


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
                        {"label": name[6:], "value": name}
                        for name in emb.columns
                        if "annot." in name
                    ]
                    + [{"label": "None", "value": "None"}],
                    value="None",
                ),
                html.Button('Clear all selections',
                    id="clear-selection",
                    n_clicks=0,
                ),
                html.Button('Undo last selection',
                    id="undo-last-selection",
                    n_clicks=0,
                ),
            ],
            style=dict(width="69%", display="inline-block"),
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
                    html.Label(html.B("Chromosome:")),
                    dcc.Dropdown(id="chrom-selector", value=chroms[0],
                                 options=options,
                                 disabled=True,
                                 style=dict(width="70%", display="inline-block"),),
                ],style={'display': 'none'}),
                html.Div([
                    html.Label(html.B("Position:")),
                    dcc.RangeSlider(
                        id="locus-selector",
                        min=0,
                        max=chromlens[chroms[0]],
                        step=1000,
                        disabled=True,
                        value=[1, 1000000],),
                ],style={'display': 'none'}),
                html.Div([
                    html.Label(html.B("Highlight:")),
                    dcc.RangeSlider(
                        id="highlight-selector", min=0, max=100, value=[25, 75],
                        disabled=True,
                    ),
                ], style={'display':'none'}),
                html.Div([
                    html.Label(html.B("Plot-type:")),
                    dcc.RadioItems(
                        id="plot-type-selector",
                        options=[
                            {"label": "line", "value": "line"},
                            {"label": "bar", "value": "bar"},
                        ],
                        value="line",
                        style=dict(width="49%", display="inline-block"),
                    ),
                ]),
                html.Div([
                    html.Label(html.B("Normalize:")),
                    dcc.RadioItems(
                        id="normalize-selector",
                        options=[
                            {"label": "No", "value": "No"},
                            {"label": "Yes", "value": "Yes"},
                        ],
                        value="No",
                        style=dict(width="49%", display="inline-block"),
                    ),
                ]),
                html.Div([
                    html.Label(html.B("Overlay tracks:")),
                    dcc.RadioItems(
                        id="overlay-track-selector",
                        options=[
                            {"label": "No", "value": "No"},
                            {"label": "Yes", "value": "Yes"},
                        ],
                        value="No",
                        style=dict(width="49%", display="inline-block"),
                    ),
                ]),
                html.Div(id="test-field", style={"display": "none"}),
            ],
            style=dict(width="29%", display="inline-block"),
        ),
        html.Br(),
        html.Div(
            [dcc.Graph(id="scatter-plot", style={'height': '500px'})],
            style=dict(width="49%", display="inline-block"),
        ),
        html.Div(
            [
                dcc.Graph(id="summary-plot", style={'height': '500px'}),
            ],
            style=dict(width="49%", display="inline-block"),
        ),
        html.Br(),
        html.Div([
                html.Div(id="stats-field"),
                dcc.Store(
                    id="session-id", data=session_id,
                ),
                dcc.Store(
                    id="dragmode-track", data="select",
                ),
                dcc.Store(
                    id="dragmode-scatter", data="select",
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



@app.callback(
    Output(component_id="chrom-selector", component_property="value"),
    [Input(component_id="gene-selector", component_property="value")],
)
@log_layer
def update_chrom_selector_value(coord):
    if coord is None:
        raise PreventUpdate
    chrom = coord.split(':')[0]
    return chrom


@app.callback(
    Output(component_id="locus-selector", component_property="max"),
    [Input(component_id="chrom-selector", component_property="value")],
)
@log_layer
def update_locus_selector_maxrange(locus):
    if locus is None:
        raise PreventUpdate
    chrom = locus.split(':')[0]
    return chromlens[chrom]


@app.callback(
    Output(component_id="locus-selector", component_property="value"),
    [Input(component_id="gene-selector", component_property="value"),
     Input(component_id="summary-plot", component_property="relayoutData"),],
)
@log_layer
def update_locus_selector_value(coord, relayout):

    ctx = dash.callback_context
    if ctx.triggered is None:
        raise PreventUpdate

    elif 'summary-plot.relayoutData' == ctx.triggered[0]['prop_id'] and \
       ctx.triggered[0]['value'] is not None and \
       'xaxis.range[0]' in ctx.triggered[0]['value']:
        interval = [
            max(0, int(relayout["xaxis.range[0]"])),
            int(relayout["xaxis.range[1]"]),
        ]
        interval.sort()
        return interval

    elif 'gene-selector.value' == ctx.triggered[0]['prop_id']:

        start, end = coord.split(':')[1].split('-')
        return [int(start), int(end)]

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
def update_scatter(annot, selection_store, dragmode, session_id):
    co = {}
    if annot in emb.columns:
        tracknames = sorted(emb[annot].unique().tolist())
    else:
        tracknames = ['None']

    colors = {trackname: annotationcolors[i%len(annotationcolors)] for i, \
              trackname in enumerate(tracknames)}

    if selection_store is not None and len(selection_store) > 0:
       selection_hash = selection_store[-1]
       selection = get_selection(session_id, selection_hash)
       selnames = [selkey for selkey in selection]
       colors.update({trackname: selectioncolors[i%len(selectioncolors)] for i, \
                 trackname in enumerate(selnames)})
       tracknames += selnames

    fig = px.scatter(
        emb,
        x="dim1",
        y="dim2",
        hover_name="barcode",
        opacity=0.3,
        color=annot if annot != "None" else None,
        color_discrete_sequence=[annotationcolors[0]] if annot == "None" else None,
        color_discrete_map=colors,
        custom_data=["barcode"],
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
    Output(component_id="summary-plot", component_property="figure"),
    [
        Input(component_id="chrom-selector", component_property="value"),
        Input(component_id="locus-selector", component_property="value"),
        Input(component_id="plot-type-selector", component_property="value"),
        Input(component_id="highlight-selector", component_property="value"),
        Input(component_id="normalize-selector", component_property="value"),
        Input(component_id="overlay-track-selector", component_property="value"),
        Input(component_id="annotation-selector", component_property="value"),
        Input(component_id="selection-store", component_property="data"),
    ],
    [
        State(component_id="dragmode-track", component_property="data"),
        State(component_id="session-id", component_property="data"),
    ],
)
@log_layer
def update_summary(chrom, interval, plottype,
                   highlight, normalize, overlay,
                   annotation, selectionstore, dragmode, session_id):
    normalize = True if normalize == "Yes" else False
    overlay = True if overlay == "No" else False

    if chrom is None:
        raise PreventUpdate
    if interval is None:
        raise PreventUpdate

    regs_ = get_region_selection(chrom, interval, highlight)

    tracknames = ['total_coverage']
    if normalize:
        regs_.loc[:,'total_coverage'] /= readsincell.sum()*1e5

    if annotation != 'None':
        df = pd.merge(cmat.cannot, emb, how='left', left_on='cell', right_on='barcode')
        tracknames = sorted(df[annotation].unique().tolist())

        for trackname in tracknames:
            cell_ids =  df[df[annotation] == trackname].index
            regs_.loc[:,trackname] = cell_coverage(cmat, regs_.index, cell_ids)
            if normalize:
                regs_.loc[:,trackname] /= readsincell[cell_ids].sum()*1e5

    colors = {trackname: annotationcolors[i%len(annotationcolors)] for i, \
              trackname in enumerate(tracknames)}

    if selectionstore is not None:
        # get cells from the cell selection store
        selnames = []
        selcells = get_cells(session_id, selectionstore[-1])
        for selcell in selcells or []:
            cell_ids = selcells[selcell]
            regs_.loc[:,selcell] = cell_coverage(cmat, regs_.index, cell_ids)
            selnames.append(selcell)
            if normalize:
                regs_.loc[:,selcell] /= readsincell[cell_ids].sum()*1e5
        colors.update({trackname: selectioncolors[i%len(selectioncolors)] for i, \
                  trackname in enumerate(selnames)})
        tracknames += selnames


    controlprops = {'plottype': plottype,
                    'dragmode': dragmode,
                    'overlay': overlay}

    trackmanager = TrackManager(regs_, tracknames,
                                colors,
                                genes, controlprops)
    fig = trackmanager.draw()

    return fig


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
                cmat.cannot[cmat.cannot.cell.isin(sel)].index.tolist()}

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
    Output(component_id="stats-field", component_property="children"),
    [
        Input(component_id="chrom-selector", component_property="value"),
        Input(component_id="locus-selector", component_property="value"),
        Input(component_id="highlight-selector", component_property="value"),
        Input(component_id="annotation-selector", component_property="value"),
        Input(component_id="selection-store", component_property="data"),
    ],
    [
     State(component_id="session-id", component_property="data"),
    ],
)
@log_layer
def update_statistics(chrom, interval,
                      highlight, annotation,
                      selectionstore, session_id):
    if chrom is None:
        raise PreventUpdate
    if interval is None:
        raise PreventUpdate

    start, end = interval
    regs_ = get_region_selection(chrom, interval, highlight)

    regs_ = regs_[regs_.highlight == "inside"]

    region_ids = regs_.index

    tablemanager = TableManager(regs_)

    if annotation != 'None':
        df = pd.merge(cmat.cannot, emb, how='left', left_on='cell', right_on='barcode')
        tracknames = sorted(df[annotation].unique().tolist())

        for trackname in tracknames:
            cell_ids =  df[df[annotation] == trackname].index
            n11 = cmat.cmat[region_ids, :][:, cell_ids].sum()
            # all selected cell reads
            c1 = readsincell[cell_ids].sum()
            r1 = readsinregion[region_ids].sum()
            n = readsincell.sum()
            tablemanager.add_row(trackname, n11, c1, r1, n, len(cell_ids))

    if selectionstore is not None:
        # get cells from the cell selection store
        selcells = get_cells(session_id, selectionstore[-1])
        for selcell in selcells or []:
            cell_ids = selcells[selcell]

            # selected cells and highlighted region
            n11 = cmat.cmat[region_ids, :][:, cell_ids].sum()
            # all selected cell reads
            c1 = readsincell[cell_ids].sum()
            r1 = readsinregion[region_ids].sum()
            n = readsincell.sum()
            tablemanager.add_row(selcell, n11, c1, r1, n, len(cell_ids))

    tab = tablemanager.draw()
    return tab


@app.callback(
    Output(component_id="highlight-selector", component_property="value"),
    [
        Input(component_id="locus-selector", component_property="value"),
        Input(component_id="summary-plot", component_property="selectedData"),
    ],
)
@log_layer
def update_highlight_selector(interval, selected):
    ctx = dash.callback_context

    if ctx.triggered is None:
        raise PreventUpdate

    if ctx.triggered[0]['prop_id'] == '.':
        #inital value
        return [25, 75]

    if interval is None or selected is None:
        raise PreventUpdate

    windowsize = interval[1] - interval[0]

    if 'range' not in selected and 'lassoPoints' not in selected:
        raise PreventUpdate
    nums = None
    if 'range' in selected:
        for key in selected['range']:
            if 'x' in key:
                nums = [int(point) for point in selected['range'][key]]
    if 'lassoPoints' in selected:
        for key in selected['lassoPoints']:
            if 'x' in key:
                nums = [int(point) for point in selected['lassoPoints'][key]]

    if nums is None:
        raise PreventUpdate
    ranges = [
        min(100, max(0, int((min(nums) - interval[0]) / windowsize * 100))),
        min(100, max(0, int((max(nums) - interval[0]) / windowsize * 100))),
    ]
    return ranges


@app.callback(
    Output(component_id="dragmode-track", component_property="data"),
    [Input(component_id="summary-plot", component_property="relayoutData"),
     ],
)
@log_layer
def update_dragmode_selector(relayout):
    if relayout is None:
        raise PreventUpdate
    if "dragmode" not in relayout:
        raise PreventUpdate
    return relayout["dragmode"]

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
