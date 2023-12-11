'''
Created on 1 mars 2023

@author: michel
'''
import os
import uuid
import tarfile
import time

import shutil
from astropy.io import fits
from constants import PATHTO

class Session(object):
    '''
    Some utilities helping to manage the session data
    '''
    
    def __init__(self, filepath):
        self.obsid=None
        self.filepath = filepath
        self.name = f"{uuid.uuid4()}"
        self.path = os.path.realpath(os.path.join(PATHTO.sessions, self.name))
        self.remove_session_directory()
        print(f"build session folder {self.path}")
        os.mkdir(self.path )
        self._set_obsid()
        
        
    def _set_obsid(self):
        """it is important to set the OBSid within the session make sure
        the is no overlap with other request since process_one_observation is not
        thread safe at all
        """
        raw_data = fits.open(self.filepath, memmap=True)
    
        obs_information = raw_data[0].header
        self.obsid = str(obs_information['OBS_ID'])
        
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
        from rest_api.logic import process_one_observation
        process_one_observation(self, self.filepath,  None)

    @staticmethod
    def clean_up(delay_hours):
        now = time.time()
        old = now - (delay_hours*3600)
        # remove older sessions
        for root, dirs, _ in os.walk(PATHTO.sessions, topdown=False):
            for _dir in dirs:
                dir_path = os.path.join(root, _dir)
                if os.path.getmtime(dir_path) < old:
                    print(f"rm session {_dir}")
                    shutil.rmtree(dir_path)     
        
        # remove older sessions
        list_precomputed_obsids = os.listdir(PATHTO.precomputed_obsids)
        if len(list_precomputed_obsids) >= 30:
            oldest_file = min(list_precomputed_obsids, key=os.path.getctime)
            print(f"rm cached file {oldest_file}")
            os.remove(os.path.abspath(oldest_file))
       
    