rm -rf cache-directory
export SCBROWSE_LOGDIR=scbrowse-logs
mkdir -p $SCBROWSE_LOGDIR

export SCBROWSE_MATRIX=data/scbrowse.h5ad
export SCBROWSE_GENES=data/danrer11_refgenes.bed
export SCBROWSE_LOGS=$SCBROWSE_LOGDIR/scbrowse.log

# debugging with:
 scbrowse

# production with:
#gunicorn -c gunicorn.config.py cli:server

