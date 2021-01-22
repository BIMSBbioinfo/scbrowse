bind = ':8000'
workers = 2
chdir = './src/scbrowse'
accesslog = '/home/sa-l-scbrowse/logs/access_gunicorn.log'
errorlog = '/home/sa-l-scbrowse/logs/error_gunicorn.log'
graceful_timeout = 60
timeout=60
#daemon=True
