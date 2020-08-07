import os

bind = ':8000'
workers = 4
chdir = './src/scbrowse'
#accesslog = '/mnt/storage/wolfgang/wolfgang/src/scbrowse-logs/access_gunicorn.log'
accesslog = os.path.join(os.environ['SCBROWSE_LOGDIR'], 'access_gunicorn.log')
errorlog = os.path.join(os.environ['SCBROWSE_LOGDIR'], 'error_gunicorn.log')
#proc_name = 'singlecell'
