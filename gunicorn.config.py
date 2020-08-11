bind = ':8000'
workers = 4
chdir = './src/scbrowse'
accesslog = '/mnt/storage/wolfgang/wolfgang/src/scbrowse-logs/access_gunicorn.log'
errorlog = '/mnt/storage/wolfgang/wolfgang/src/scbrowse-logs/error_gunicorn.log'
graceful_timeout = 60
timeout=60
daemon = True
