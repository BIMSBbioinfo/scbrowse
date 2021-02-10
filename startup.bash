if [ -f "/home/sa-l-scbrowse/miniconda3/etc/profile.d/conda.sh" ]; then
  source /home/sa-l-scbrowse/miniconda3/etc/profile.d/conda.sh
fi
source /home/sa-l-scbrowse/miniconda3/etc/profile.d/conda.sh

/home/sa-l-scbrowse/miniconda3/bin/conda init bash
/home/sa-l-scbrowse/miniconda3/bin/conda activate scbrowse

rm -rf cache-directory
export SCBROWSE_LOGDIR=/home/sa-l-scbrowse/logs
mkdir -p $SCBROWSE_LOGDIR

export SCBROWSE_MATRIX=/home/sa-l-scbrowse/data/scbrowse.h5ad
export SCBROWSE_GENES=/home/sa-l-scbrowse/data/danrer11_refgenes.bed
export SCBROWSE_LOGS=$SCBROWSE_LOGDIR/scbrowse.log

# debugging with:
#scbrowse

# production with:
/home/sa-l-scbrowse/miniconda3/envs/scbrowse/bin/gunicorn -c gunicorn.config.py cli:server

