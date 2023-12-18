'''
Created on 1 mars 2023

@author: michel
'''
import os
import uuid
import tarfile
import time
import traceback
import sys
import shutil
from multiprocessing import Process, Queue
from astropy.io import fits
from flask import send_from_directory
from werkzeug.utils import secure_filename
from constants import PATHTO

class Session(object):
    '''
    Some utilities helping to manage the session data
    '''
    
    def __init__(self, file):
        self.obsid=None
        self.filename = secure_filename(file.name)
        self.name = f"{uuid.uuid4()}"
        self.path = os.path.realpath(os.path.join(PATHTO.sessions, self.name))
        self.remove_session_directory()
        print(f"build session folder {self.path}")
        os.mkdir(self.path )
        self.filepath  = os.path.join(self.path, self.filename)
        print(f"save {self.filename} in {self.filepath}")
        file.save(self.filepath)
        self._set_obsid()
        
        
    def _set_obsid(self):
        """it is important to set the OBSid within the session make sure
        the is no overlap with other request since process_one_observation is not
        thread safe at all
        """
        raw_data = fits.open(self.filepath, memmap=True)  
        obs_information = raw_data[0].header
        self.obsid = str(obs_information['OBS_ID'])
        raw_data.close()
        
        
    def remove_session_directory(self):
        for root, dirs, _ in os.walk(PATHTO.sessions):
            for d in dirs:
                if d == self.name:
                    print(f"Remove session folder {self.name}")
                    shutil.rmtree(os.path.join(root, d))     
     
    def get_tarball(self):  
        if self.obsid is not None: 
            tarfile_prefix = self.obsid
        else:
            tarfile_prefix = "alerts"

        pdfs = []
        for file in os.listdir(self.path):
            if file.endswith(".PDF"):
                pdfs.append(file)
   
        output = os.path.join(self.path, f"{tarfile_prefix}.tgz")
        print(f"packing the result in {output}")
        with tarfile.open(output, "w:gz") as tar:
            for pdf in pdfs:
                tar.add(os.path.join(self.path, pdf), arcname=pdf)
        return output

    def store_source_list(self, obsmli_path):
        shutil.copy(obsmli_path, self.path)

    def process_observation(self):
        """run the processinf
        local import to avid circular imports Session<->logic
        """
        try:
            from rest_api.logic import process_one_observation
            q = Queue()
            p = Process(target=process_one_observation, args=(self, q))
            p.start()
            result = q.get()
            p.join()
            if result["status"] == "failed":
                return result, 500
            
            if result["nb_alerts"] == "0":
                return result, 404
            tarball_path = self.get_tarball()
            directory, filename = os.path.split(tarball_path)
            return send_from_directory(directory, filename)  , 200      

        except Exception as exp:
            traceback.print_exc(file=sys.stdout)
            result = {'status': 'failed',
              'message': f"Something went wrong {exp}"}
            return result, 500
       
       
        result = {'status': 'failed',
                  'message': f"Prohibited filename {self.filename}"}
        return result, 500

        
      
    