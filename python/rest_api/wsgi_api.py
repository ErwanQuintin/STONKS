"""
Created on March 2023

Launcher for the REST service in production  mode (wsgi)
Must be run with gunicorn: ``gunicorn --bind 0.0.0.0:5000 api.wsgi_api:application``

@author: michel
"""
import os
from rest_api.flask_app import app as application
from cache import Cache

print("starting cache monitoring")
Cache.running = True
Cache().start()    
print("cache monitoring started")
print("Mplt backend ",  os.environ["MPLBACKEND"])
if __name__ == "__main__":
    application.run()
