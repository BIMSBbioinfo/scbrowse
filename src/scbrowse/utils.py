import time
import base64
import io

import numpy as np
from coolbox.api import *
from coolbox.utilities import split_genome_range
import scanpy as sc
import matplotlib.pyplot as plt

class SingleCellTrack(HistBase):

    @classmethod
    def totalbulk(cls, frame, adata, **kwargs):
        f = cls(adata.obs, 'total', title='total',
                color='black',
                #max_value=max_value,
                **kwargs) + TrackHeight(.5)
        if frame is None:
            frame = f
        else:
            frame += f
        return frame

    @classmethod
    def pseudobulk(cls, frame, adata, groupby,
                   palettes=sc.pl.palettes.vega_20_scanpy,
                   **kwargs):
        if groupby in adata.var.columns:
            names = sorted(adata.var[groupby].unique().tolist())

        colorname = groupby + '_colors'
        if colorname not in adata.uns:
            adata.uns[colorname] = palettes[:len(names)]

        obs = adata.obs
        max_value = obs.loc[:,[f'{groupby}_{name}' for name in names]].values.max()

        for i, name in enumerate(names):

            f = cls(obs, f'{groupby}_{name}', title=f'{name}',
                    color=adata.uns[colorname][i],
                    max_value=max_value,
                    **kwargs) + TrackHeight(.5)
            if frame is None:
                frame = f
            else:
                frame += f
            frame += HLine()
        return frame


    @classmethod
    def selection(cls, frame, adata, selections,
                   palettes=sc.pl.palettes.zeileis_28,
                   **kwargs):
        names = [k for k in selections]

        colorname = 'selection_colors'
        if colorname not in adata.uns:
            adata.uns[colorname] = palettes[:len(names)]

        obs = adata.obs.copy()
        for i, name in enumerate(selections):
            sadata = adata[:,adata.var.index.isin(selections[name])]
            da = np.asarray(sadata.X.sum(1)).flatten() *1e5 / sadata.var.nFrags.sum()
            obs.loc[:,name] = da

        max_value = obs.loc[:,[name for name in selections]].values.max()
        for i, name in enumerate(selections):
            f = cls(obs, name, title=f'{name}',
                    color=adata.uns[colorname][i],
                    max_value=max_value,
                    **kwargs) + TrackHeight(.9)
            if frame is None:
                frame = f
            else:
                frame += f
            frame += HLine()
        return frame


    def __init__(self, obs, group, **kwargs):
        properties = HistBase.DEFAULT_PROPERTIES.copy()
        properties.update({
            'type': properties['style'],
            "file": 'dummy',
            #'style': 'fill',
            **kwargs,
        })
        super().__init__(**properties)
        self.obs = obs
        self.group = group

    def fetch_data(self, gr, **kwargs):
        chrom, start, end = split_genome_range(gr)
        obs = self.obs
        obs = obs.loc[(obs.chrom==chrom) & (obs.start>=start) & (obs.end<=end),self.group]
        return obs.values




def draw_figure(inputs):
    adata, locus, annotation, selcells, genefile, lock = inputs

    st=time.time()
    chrom, start, end = split_genome_range(locus)
    frame = Frame(width=15, title=chrom, fontsize=5)
    frame += XAxis(fontsize=7, title=chrom)
    sada = adata[(adata.obs.chrom==chrom) & (adata.obs.start>=start) & (adata.obs.end<=end),:]
    frame = SingleCellTrack.totalbulk(frame, sada, fontsize=5)

    if annotation != "None":
        frame = SingleCellTrack.pseudobulk(frame, sada, annotation, fontsize=5)

    if selcells is not None:
        frame = SingleCellTrack.selection(frame, sada, selcells, fontsize=5)

    frame += Spacer(.3)
    frame += BED(genefile, title='Genes')
    print(f'  finished preparation: {time.time() - st} ')
    lock.acquire()
    st=time.time()
    fig = frame.plot(f'{chrom}:{start}-{end}')
    print(f'  finished plot: {time.time() - st} ')
    st=time.time()
    buf = io.BytesIO()
    fig.savefig(buf, format = "png", bbox_inches = "tight")
    plt.close(fig)
    lock.release()
    print(f'  finished savefig: {time.time() - st} ')
    data = base64.b64encode(buf.getbuffer()).decode("ascii")
    return f"data:image/png;base64,{data}"

