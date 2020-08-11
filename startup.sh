rm -rf cache-directory
export SCBROWSE_LOGDIR=/home/sa-l-scbrowse/logs
mkdir -p $SCBROWSE_LOGDIR

export SCBROWSE_EMBEDDING=/home/sa-l-scbrowse/data/cellannot.tsv
export SCBROWSE_MATRIX=/home/sa-l-scbrowse/data/all.npz
export SCBROWSE_REGIONS=/home/sa-l-scbrowse/data/tile.bed
export SCBROWSE_GENES=/home/sa-l-scbrowse/data/danrer11_refgenes.bed
export SCBROWSE_LOGS=$SCBROWSE_LOGDIR/scbrowse.log

# debugging with:
#scbrowse

# production with:
gunicorn -c gunicorn.config.py cli:server

