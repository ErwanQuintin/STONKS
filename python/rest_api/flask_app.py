'''
Created on 1 mars 2023

@author: michel
'''
import os
from flask import Flask, request, send_file
from constants import PATHTO
from session import Session
from cache import Cache

UPLOAD_FOLDER = '/path/to/the/uploads'
ALLOWED_EXTENSIONS = ['zip','ZIP','tar','TAR','gz','GZ', 'fit', 'fits', 'FITS', 'FIT']
app = Flask(__name__)
app.secret_key = "super secret key"
app.config['MAX_CONTENT_LENGTH'] = 1024 * 1024    # 1 Mb limit
  
def allowed_file(filename):    
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.route('/stonks/process-obs', methods=['POST'])
def upload_file():
    if request.method == 'POST':

        # check if the post request has the file part
        if 'file' not in request.files:
            result = {'status': 'ko',
                      'message': 'No file part in the request'}
            return result, 400
        file = request.files['file']
        # If the user does not select a file, the browser submits an
        # empty file without a filename.
        if file.filename == '':
            result = {'status': 'ko',
                      'message': 'No selected file'}
            return result, 400
        Session.clean_up()
        return Session(file).process_observation()
                   
@app.route("/stonks/doc")
def doc():

    return send_file(os.path.join(PATHTO.basedir, 'doc', 'STONKS_Documentation.pdf'))

@app.route('/stonks//static/<path:subpath>')
def get_static(subpath):
    return send_file(os.path.join(PATHTO.basedir, 'static', subpath))

@app.route('/stonks/doc/<path:subpath>')
def get_doc(subpath):
    return send_file(os.path.join(PATHTO.basedir, 'doc', subpath))

@app.route("/stonks/")
def home():
     
    return send_file(os.path.join(PATHTO.basedir, 'doc', "index2.html"))


if __name__ == "__main__":
    """
    Test command:  curl -J -X POST -F file=@P0804670301EPX000OBSMLI0000.FIT http://127.0.0.1:5000/process-obs
    """
    print("starting cache monitoring")
    Cache.running = True
    Cache().start()
    app.run(host="127.0.0.1", port=5000)

