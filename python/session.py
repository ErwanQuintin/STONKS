'''
Created on 1 mars 2023

@author: michel
'''
import os
import uuid
import tarfile
import time

import shutil
from constants import PATHTO

class Session(object):
    '''
    Some utilities helping to manage the session data
    '''
    
    def __init__(self):
        self.obsid=None
        self.name = f"{uuid.uuid4()}"
        self.path = os.path.realpath(os.path.join(PATHTO.sessions, self.name))
        self.remove_session_directory()
        print(f"build session folder {self.path}")
        os.mkdir(self.path )
        
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
        
    