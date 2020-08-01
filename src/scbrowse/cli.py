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
import argparse
from functools import wraps
import logging
#from copy import copy
import argparse
import dash  # pylint: disable=import-error
from dash.exceptions import PreventUpdate
import dash_core_components as dcc  # pylint: disable=import-error
import dash_html_components as html  # pylint: disable=import-error
import plotly.graph_objs as go  # pylint: disable=import-error
import plotly.express as px
from plotly.subplots import make_subplots
from scipy.stats import fisher_exact

from scipy.sparse import diags
import numpy as np
from dash.dependencies import Input, Output, State
from pybedtools import BedTool, Interval

import pandas as pd
from scregseg.countmatrix import CountMatrix


parser = argparse.ArgumentParser(description="single-cell Explorer")

parser.add_argument(
    "-embedding", dest="embedding", type=str,
    help="Table with 2D embedding of cells"
)
parser.add_argument("-matrix", dest="matrix", type=str, help="Count matrix")
parser.add_argument("-bed", dest="bed", type=str, help="BED file")
parser.add_argument(
    "-genes", dest="genes", type=str, help="Gene annotation in BED12 format"
)
parser.add_argument(
    '-port', dest='port', type=int,
    default=8051,
    help="Port to use for the web app"
)
parser.add_argument(
    '-log', dest='log', type=str,
    default='scbrower.log',
    help="Port to use for the web app"
)

#args = parser.parse_args()

#parser = argparse.ArgumentParser(description='Command description.')
#parser.add_argument('names', metavar='NAME', nargs=argparse.ZERO_OR_MORE,
#                    help="A name of something.")

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
        m = m[:, cmat.cannot[cmat.cannot.cell.isin(cells)].index]
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
    #sregs_ = regs_[regs_.highlight == "inside"]
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

class TrackManager:
    def __init__(self, regs_, tracknames,
                 genes, controlprops):
        self.regs_ = regs_
        self.tracknames = tracknames
        self.trackheight = 3
        self.locus = get_locus(regs_)
        self.highlight = get_highlight(regs_)
        self.controlprops = controlprops
        #self.plottype = plottype
        #self.separated = separated
        self.genes = genes
        self.init_trace()

    def allocate(self):
        if not self.controlprops['separated']:
            ntracks = 1
        else:
            ntracks = len(self.tracknames)
       #     print('not separated')

       #     specs = [[{"rowspan": self.trackheight}]] + [[None]] * (self.trackheight-1) + [[{"rowspan": 1}]]

       #     fig = make_subplots(
       #         rows=ntracks * self.trackheight + 1,
       #         cols=1,
       #         specs=specs,
       #         shared_xaxes=True,
       #         vertical_spacing=0.05,
       #         print_grid=True,
       #     )
       # else:
        #print(' separated')
        #ntracks = len(self.tracknames)
        specs=[]
        for n in range(ntracks):
            specs += [[{"rowspan": self.trackheight}]]
            specs += [[None]]*(self.trackheight - 1)
        specs += [[{"rowspan": 1}]]
        #print(specs)
        #print(len(specs))
        #specs = [[{"rowspan": ntracks * 3}]] + \
        #    [[None]] * (ntracks * 3 + 1) + [[{"rowspan": 1}]]

        fig = make_subplots(
            rows=ntracks * self.trackheight + 1,
            cols=1,
            specs=specs,
            shared_xaxes=True,
            vertical_spacing=0.05,
            print_grid=True,
        )
        return fig

    def draw_summary_track(self, trackname):
        plobjs = []
        plottype = self.controlprops['plottype']
        sregs_ = self.regs_
        if plottype == "bar":
            plobjs.append(go.Bar(x=sregs_.start,
                                 y=sregs_[trackname],
                                 name=trackname))
        else:
            plobjs.append(go.Scatter(x=sregs_.start,
                                     y=sregs_[trackname],
                                     mode="lines", name=trackname,))

        plobjs.append(go.Scatter(
                x=sregs_.start,
                y=sregs_[trackname],
                mode="markers",
                opacity=0,
                hoverinfo="skip",
                showlegend=False,
            ))
        return plobjs

    def init_trace(self):
        self.itrace = 1

    def next_trace(self):
        if self.controlprops['separated']:
            self.itrace += self.trackheight

    def extend_trace(self, fig, plobjs):
        for plobj in plobjs:
            print(f'extend_trace at {self.itrace}')
            fig.add_trace(plobj, row=self.itrace, col=1)
        return fig

    def draw_tracks(self, fig):
        self.init_trace()
        for trackname in self.tracknames:
            self.extend_trace(fig, self.draw_summary_track(trackname))
            self.next_trace()

        #hl = get_highlight(self.regs_)

        #if selected is None:
        #    fig = draw_summary_figure_3(fig, regs_, plottype, "total_coverage")
        #else:
        #    fig = draw_summary_figure_3(fig, regs_, plottype, "selected")
        #    fig = draw_summary_figure_3(fig, regs_, plottype, "total_coverage")

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
        #return self.draw_highlight(fig, regs_)

    def draw_annotation(self, fig):
        if self.genes is not None:
            ntracks = len(self.tracknames) if self.controlprops['separated'] else 1
            chrom, start, end = self.locus
            plobjs = _draw_gene_annotation(fig, self.genes, chrom, start, end)
            for plobj in plobjs:
                fig.add_trace(
                    plobj,
                    row=ntracks*self.trackheight + 1,
                    col=1,
                )
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

def draw_highlight(fig, regs_):
    hl = get_highlight(regs_)
    if hl is None:
        return fig
    xmin, xmax = hl

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


def draw_summary_figure_3(fig, regs_, plottype, yfield):
    if plottype == "bar":
        fig.add_trace(
            go.Bar(x=regs_.start, y=regs_[yfield], name=yfield), row=1, col=1,
        )
    else:
        fig.add_trace(
            go.Scatter(x=regs_.start, y=regs_[yfield], mode="lines", name=yfield,),
            row=1,
            col=1,
        )

    fig.add_trace(
        go.Scatter(
            x=regs_.start,
            y=regs_[yfield],
            mode="markers",
            opacity=0,
            hoverinfo="skip",
            showlegend=False,
        ),
        row=1,
        col=1,
    )
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
#    else:
#        return []
#        fig.add_trace(
#            go.Scatter(
#                x=xs,
#                y=ys,
#                mode="lines",
#                fill="toself",
#                name="Genes",
#                marker=dict(color="goldenrod"),
#            ),
#            row=12,
#            col=1,
#        )
#        fig.add_trace(
#            go.Scatter(
#                x=midpoints,
#                y=[0.5] * len(midpoints),
#                text=names,
#                mode="text",
#                opacity=0.0,
#                name="Genes",
#                customdata=rangeannot,
#                hovertemplate="%{text}<br>%{customdata}",
#                showlegend=False,
#            ),
#            row=12,
#            col=1,
#        )
#        fig.layout["yaxis2"]["showticklabels"] = False
#        fig.layout["yaxis2"]["showgrid"] = False
#        fig.layout["yaxis2"]["zeroline"] = False
#        fig.layout["xaxis2"]["showgrid"] = False
#        fig.layout["xaxis2"]["zeroline"] = False
#    else:
#        fig.layout["xaxis"]["showticklabels"] = True
#    return fig

def draw_gene_annotation(fig, genes, chrom, start, end):
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
    #    plobjs = [
    #        go.Scatter(
    #            x=xs,
    #            y=ys,
    #            mode="lines",
    #            fill="toself",
    #            name="Genes",
    #            marker=dict(color="goldenrod"),
    #        ),
    #        go.Scatter(
    #            x=midpoints,
    #            y=[0.5] * len(midpoints),
    #            text=names,
    #            mode="text",
    #            opacity=0.0,
    #            name="Genes",
    #            customdata=rangeannot,
    #            hovertemplate="%{text}<br>%{customdata}",
    #            showlegend=False,
    #        ),
    #    ]
    #    return plobjs
    #else:
    #    return []
        fig.add_trace(
            go.Scatter(
                x=xs,
                y=ys,
                mode="lines",
                fill="toself",
                name="Genes",
                marker=dict(color="goldenrod"),
            ),
            row=12,
            col=1,
        )
        fig.add_trace(
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
            row=12,
            col=1,
        )
        fig.layout["yaxis2"]["showticklabels"] = False
        fig.layout["yaxis2"]["showgrid"] = False
        fig.layout["yaxis2"]["zeroline"] = False
        fig.layout["xaxis2"]["showgrid"] = False
        fig.layout["xaxis2"]["zeroline"] = False
    else:
        fig.layout["xaxis"]["showticklabels"] = True
    return fig



def main(args=None):
    args = parser.parse_args(args=args)
    #args = parser.parse_args()

    logging.basicConfig(filename = args.log,
                        level=logging.DEBUG,
                        format='%(asctime)s;%(levelname)s;%(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')

    logging.debug('scbrowser - startup')
    logging.debug(args)

    ##############
    # load data
    ##############

    genes = BedTool(args.genes)

    logging.debug(f'Number of genes: {len(genes)}')

    emb = pd.read_csv(args.embedding, sep="\t")
    if "barcode" not in emb.columns:
        emb["barcode"] = emb.index

    logging.debug(f'Embedding dims: {emb.shape}')

    cmat = CountMatrix.from_mtx(args.matrix, args.bed)

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

    genelocus = [dict(label=g.name, value=f'{g.chrom}:{g.start}-{g.end}') for g in genes]

    ##############
    # instantiate app
    ##############

    app = dash.Dash("scbrowser")

    app.title = "scbrowser"
    app.config["suppress_callback_exceptions"] = True

    ##############
    # define layout
    ##############

    app.layout = html.Div(
        [
            html.Div(
                [
                    html.Div(
                        id="dragmode-selector", children="select", style={"display": "none"}
                    ),
                    html.Label(html.B("Annotation:")),
                    dcc.Dropdown(
                        id="annot-selector",
                        options=[
                            {"label": name[6:], "value": name}
                            for name in emb.columns
                            if "annot." in name
                        ]
                        + [{"label": "None", "value": "None"}],
                        value="None",
                    ),
                    html.Div(id="stats-field"),
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
                        html.Label(html.B("Chromosome:")),
                        dcc.Dropdown(id="chrom-selector", value=chroms[0],
                                     options=options,
                                     style=dict(width="70%", display="inline-block"),),
                    ],),
                    html.Div([
                        html.Label(html.B("Position:")),
                        dcc.RangeSlider(
                            id="locus-selector",
                            min=0,
                            max=chromlens[chroms[0]],
                            step=1000,
                            value=[1, 1000000],),
                    ],),
                    html.Div([
                        html.Label(html.B("Highlight:")),
                        dcc.RangeSlider(
                            id="highlight-selector", min=0, max=100, value=[25, 75],
                        ),
                    ]),
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
                        html.Label(html.B("Separate tracks:")),
                        dcc.RadioItems(
                            id="separate-track-selector",
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
                style=dict(width="49%", display="inline-block"),
            ),
            html.Br(),
            html.Div(
                [dcc.Graph(id="scatter-plot")],
                style=dict(width="49%", display="inline-block"),
            ),
            html.Div(
                [
                    dcc.Graph(id="summary-plot"),
                    # html.Pre(id="stats-field"),
                ],
                style=dict(width="49%", display="inline-block"),
            ),
        ],
        className="container",
    )


    ############
    # app callbacks
    ############



    @app.callback(
        Output(component_id="chrom-selector", component_property="value"),
        [Input(component_id="gene-selector", component_property="value")],
    )
    @log_layer
    def update_chrom_selector_value(coord):
        chrom = coord.split(':')[0]
        return chrom




    @app.callback(
        Output(component_id="locus-selector", component_property="max"),
        [Input(component_id="chrom-selector", component_property="value")],
    )
    @log_layer
    def update_locus_selector_maxrange(locus):
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
        [Input(component_id="annot-selector", component_property="value")],
    )
    @log_layer
    def update_scatter(annot):
        fig = px.scatter(
            emb,
            x="dim1",
            y="dim2",
            hover_name="barcode",
            opacity=0.3,
            color=annot if annot != "None" else None,
            custom_data=["barcode"],
            template="plotly_white",
        )

        # fig.update_data
        fig.update_layout(clickmode="event+select")

        return fig


    @app.callback(
        Output(component_id="summary-plot", component_property="figure"),
        [
            Input(component_id="chrom-selector", component_property="value"),
            Input(component_id="locus-selector", component_property="value"),
            Input(component_id="plot-type-selector", component_property="value"),
            Input(component_id="scatter-plot", component_property="selectedData"),
            Input(component_id="highlight-selector", component_property="value"),
            Input(component_id="dragmode-selector", component_property="children"),
            Input(component_id="normalize-selector", component_property="value"),
            Input(component_id="separate-track-selector", component_property="value"),
        ],
    )
    @log_layer
    def update_summary(chrom, interval, plottype, selected,
                       highlight, dragmode, normalize, separated):
        normalize = True if normalize == "Yes" else False
        separated = True if separated == "Yes" else False

        if chrom is None:
            raise PreventUpdate
        if interval is None:
            raise PreventUpdate

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

        tracknames = ['total_coverage']
        if normalize:
            regs_['total_coverage'] /= readsincell.sum()*1e5

        if selected is not None:
            sel = [point["customdata"][0] for point in selected["points"]]
            cell_ids = cmat.cannot[cmat.cannot.cell.isin(sel)].index
            regs_["selected"] = cell_coverage(cmat, regs_.index, sel)
            tracknames += ['selected']

            if normalize:
                regs_['selected'] /= readsincell[cell_ids].sum()*1e5

        controlprops = {'plottype': plottype,
                        'dragmode': dragmode,
                        'separated': separated}
        trackmanager = TrackManager(regs_, tracknames,
                                    genes, controlprops)
        fig = trackmanager.draw()

        #specs = [[{"rowspan": 11}]] + [[None]] * 10 + [[{"rowspan": 1}]]

        #fig = make_subplots(
        #    rows=12,
        #    cols=1,
        #    specs=specs,
        #    shared_xaxes=True,
        #    vertical_spacing=0.01,
        #    print_grid=True,
        #)

        #if selected is None:
        #    fig = draw_summary_figure_3(fig, regs_, plottype, "total_coverage")
        #else:
        #    fig = draw_summary_figure_3(fig, regs_, plottype, "selected")
        #    fig = draw_summary_figure_3(fig, regs_, plottype, "total_coverage")

        #hl = get_highlight(regs_)
        #title = f"Window: {chrom}:{start}-{end}"
        #if hl is not None:
        #    title += f"<br>Highlight: {chrom}:{hl[0]}-{hl[1]}"

        #fig.layout.update(
        #    dict(title=title, dragmode=dragmode, clickmode="event+select",
        #         template="plotly_white")
        #)
        ##fig.update_layout(clickmode="event+select")
        ##fig.update_layout(template="plotly_white")

        #fig = draw_highlight(fig, regs_)

        #if genes is not None:
        #    fig = draw_gene_annotation(fig, genes, chrom, start, end)


        return fig


    @app.callback(
        Output(component_id="stats-field", component_property="children"),
        [
            Input(component_id="chrom-selector", component_property="value"),
            Input(component_id="locus-selector", component_property="value"),
            Input(component_id="scatter-plot", component_property="selectedData"),
            Input(component_id="highlight-selector", component_property="value"),
        ],
    )
    @log_layer
    def update_statistics(chrom, interval, selected, highlight):
        if chrom is None:
            raise PreventUpdate
        if interval is None:
            raise PreventUpdate

        start, end = interval
        regs_ = regs.copy()

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

        regs_ = regs_[regs_.highlight == "inside"]

        region_ids = regs_.index

        if selected is not None:
            sel = [point["customdata"][0] for point in selected["points"]]
            cell_ids = cmat.cannot[cmat.cannot.cell.isin(sel)].index

            # selected cells and highlighted region
            n11 = cmat.cmat[region_ids, :][:, cell_ids].sum()
            # all selected cell reads
            c1 = readsincell[cell_ids].sum()

        else:
            n11 = 0
            c1 = 0
            cell_ids = []

        N = readsincell.sum()
        # all non-selected cell reads
        c2 = N - c1

        # all highlighted reads
        r1 = readsinregion[region_ids].sum()
        # all non-highlighted reads
        r2 = N - readsinregion[region_ids].sum()

        # selected cells and non-highlighted reads
        n21 = c1 - n11
        # non-selected cells and highlighted reads
        n12 = r1 - n11
        # non-selected cells and non-highlighted reads
        n22 = N - n11 - n12 - n21

        def get_ele(n):
            return html.Td(n, style=dict(textAlign='right',
                                         backgroundColor='lightblue'))

        def get_mele(n):
            return html.Td(n, style=dict(textAlign='right',
                                         backgroundColor='lightgray'))

        tab = html.Table(
            [
                html.Thead([html.Tr([html.Th("Cells", colSpan=3)]),]),
                html.Thead(
                    [
                        html.Tr(
                            [
                                html.Th(f"selected cells: {len(cell_ids)}"),
                                html.Th(
                                    "remaining cells: {}".format(
                                        cmat.cmat.shape[1] - len(cell_ids)
                                    )
                                ),
                                html.Th(f"total cells: {cmat.cmat.shape[1]}"),
                            ]
                        ),
                    ]
                ),

                html.Tbody(
                    [
                        html.Tr(
                            [
                                get_ele(n11),
                                get_ele(n12),
                                get_mele(r1),
                                html.Td("No. reads in highlight"),
                            ]
                        ),
                        html.Tr(
                            [
                                get_ele(n21),
                                get_ele(n22),
                                get_mele(r2),
                                html.Td("No. reads outside highlight"),
                            ]
                        ),
                        html.Tr(
                            [get_mele(c1), get_mele(c2), get_mele(N), html.Td(""),]
                        ),
                    ]
                ),
            ]
        )

        odds, pval = fisher_exact([[n11, n12], [n21, n22]], alternative='greater')
        return html.Div(
            [
                tab,
                html.Pre(
                    children=f"""Fisher's test:
                      P-value={pval}
                      odds={odds}"""
                ),
            ]
        )


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

        if 'range' in selected:
            nums = [int(point) for point in selected['range']['x']]
        if 'lassoPoints' in selected:
            nums = [int(point) for point in selected['lassoPoints']['x']]

        ranges = [
            min(100, max(0, int((min(nums) - interval[0]) / windowsize * 100))),
            min(100, max(0, int((max(nums) - interval[0]) / windowsize * 100))),
        ]
        return ranges


    @app.callback(
        Output(component_id="dragmode-selector", component_property="children"),
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

    #@app.callback(
    #    Output(component_id="test-field", component_property="children"),
    #    [
#Inp#ut(component_id="scatter-plot", component_property="relayoutData"),
    #     Input(component_id="scatter-plot", component_property="selectedData"),
    #     Input(component_id="test-field", component_property="children")
    #     ],
    #)
    #@log_layer
    #def update_point_selector(selected, data):
    #    print(selected, data)
    #    return 'hallo'


    ############
    # run server
    ############

    app.run_server(debug=True, port=args.port)
