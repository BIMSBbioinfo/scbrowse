rm -rf cache-directory
export SCBROWSE_LOGDIR=/mnt/storage/wolfgang/wolfgang/src/scbrowse-logs
mkdir -p $SCBROWSE_LOGDIR

export SCBROWSE_EMBEDDING=/mnt/storage/wolfgang/wolfgang/src/python-scbrowse/data/cellannot.tsv
export SCBROWSE_MATRIX=/mnt/storage/wolfgang/wolfgang/src/python-scbrowse/data/all.npz
export SCBROWSE_REGIONS=/mnt/storage/wolfgang/wolfgang/src/python-scbrowse/data/tile.bed
export SCBROWSE_GENES=/mnt/storage/wolfgang/wolfgang/src/python-scbrowse/data/danrer11_refgenes.bed
export SCBROWSE_LOGS=$SCBROWSE_LOGDIR/scbrowse.log

# debugging with:
# scbrowse

# production with:
gunicorn -c gunicorn.config.py cli:server

