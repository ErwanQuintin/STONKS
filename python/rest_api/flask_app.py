'''
Created on 1 mars 2023

@author: michel
'''
import os
import json
from flask import Flask, request, send_from_directory, send_file, jsonify
from werkzeug.utils import secure_filename
from multiprocessing import Process, Queue
from constants import PATHTO
from rest_api.logic import process_one_observation
from session import Session

UPLOAD_FOLDER = '/path/to/the/uploads'
ALLOWED_EXTENSIONS = ['fit', 'fits']
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
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            try:
                session = Session()
                obsmli_path = os.path.join(session.path, filename)
                file.save(obsmli_path)
                q = Queue()
                p = Process(target=session.process_observation, args=(q))
                p.start()
                result = q.get()
                p.join()
                if result["status"] == "failed":
                    return result, 500
                
                if result["nb_alerts"] == "0":
                    return result, 404
                session.obsid = result["obsid"]
                tarball_path = session.get_tarball()
                directory, filename = os.path.split(tarball_path)
                Session.clean_up(240)
                return send_from_directory(directory, filename)  , 200      

            except Exception as exp:
                result = {'status': 'failed',
                  'message': f"Something went wrong {str(exp)}"}
                return result, 500
       
       
        result = {'status': 'failed',
                  'message': f"Prohibited filename {file.filename}"}
        return result, 500

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
    app.run(host="0.0.0.0", port=5000)

