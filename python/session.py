'''
Created on 1 mars 2023

@author: michel
'''
import os
import uuid
import tarfile
import traceback
import sys
import shutil
from pathlib import Path
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
        self.obsmli_path=None
        self.image_path=None
        self.filename = secure_filename(file.filename)
        self.name = f"{uuid.uuid4()}"
        self.path = os.path.realpath(os.path.join(PATHTO.sessions, self.name))
        self.remove_session_directory()
        print(f"build session folder {self.path}")
        os.mkdir(self.path )
        self.inputtarpath  = os.path.join(self.path, self.filename)
        print(f"save {self.filename} in {self.inputtarpath}")
        file.save(self.inputtarpath)

        if self.filename.endswith('.FIT'):
            print('old version with only OBMSLI, no image found')
            self.obsmli_path = os.path.join(self.path, self.filename)
        else:
            print(f'extract the input data in {self.inputtarpath}')
            with tarfile.open(self.inputtarpath, 'r') as tar:
                tar.extractall(self.path, filter='data')
            print('initializing and checking the input file paths')
            for file in os.listdir(self.path):
                if 'OBSMLI' in file:
                    self.obsmli_path = os.path.join(self.path, file)
                elif 'IMAGE' in file:
                    self.image_path = os.path.join(self.path, file)


    def _set_obsid(self):
        """it is important to set the OBSid within the session make sure
        the is no overlap with other request since process_one_observation is not
        thread safe at all
        """
        try :
            raw_data = fits.open(self.obsmli_path, memmap=True)
            obs_information = raw_data[0].header
            self.obsid = str(obs_information['OBS_ID'])
        except:
            return False
        raw_data.close()
        return True

        
        
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
            if file.endswith(".PDF") or file.endswith('.json'):
                pdfs.append(file)
   
        output = os.path.join(self.path, f"{tarfile_prefix}.tgz")
        print(f"packing the result in {output}")
        
        with tarfile.open(output, "w:gz") as tar:
            for pdf in pdfs:
                import psutil
                print(pdf)
                print("CPU usage (%):", psutil.cpu_percent(interval=1))
        
                ram = psutil.virtual_memory()
                print("RAM usage (%):", ram.percent)
                print("RAM used (GB):", round(ram.used / 1e9, 2))        
                total, used, free = map(int, os.popen('free -t -m').readlines()[-1].split()[1:])
                print("free ", round((used / total) * 100, 2))
                tar.add(os.path.join(self.path, pdf), arcname=pdf)
                total, used, free = map(int, os.popen('free -t -m').readlines()[-1].split()[1:])
                print("free out ", round((used / total) * 100, 2))
        return output

    def store_source_list(self, obsmli_path):
        shutil.copy(obsmli_path, self.path)

    def process_observation(self):
        """run the processinf
        local import to avid circular imports Session<->logic
        """
        try:
            
            if not self._set_obsid():
                result = {'status': "failed",
                  'message': f"{self.obsmli_path} does not look like a OBSMLI fits file"}
                return result, 400

            from rest_api.logic import process_one_observation
            q = Queue()
            p = Process(target=process_one_observation, args=(self, q))
            p.start()
            result = q.get()
            p.join()
            if result["status"].startswith("failed"):
                return result, 500
            
            if result["nb_alerts"] == "0":
                return result, 404
            tarball_path = self.get_tarball()
            directory, filename = os.path.split(tarball_path)
                
            total, used, free = map(int, os.popen('free -t -m').readlines()[-1].split()[1:])
            print("free out ", total, " " , used, " " , free)

            return send_from_directory(directory, filename)  , 200      

        except Exception as exp:
            traceback.print_exc(file=sys.stdout)
            result = {'status': f"failed on obs {self.obsid}",
              'message': f"Something went wrong {exp}"}
            return result, 500
              
        result = {'status': "failed",
                  'message': f"Prohibited filename {self.filename}"}
        return result, 500

        
    @staticmethod
    def clean_up():
        folder = Path(PATHTO.sessions)
        if not folder.is_dir():
            raise ValueError(f"{PATHTO.sessions} is not a valid directory")
    
        # Get all subdirectories
        subfolders = [f for f in folder.iterdir() if f.is_dir()]
    
        # Sort by last modification time (most recent first)
        subfolders.sort(key=lambda f: f.stat().st_mtime, reverse=True)
    
        # Keep only the 20 most recent
        folders_to_delete = subfolders[40:]
    
        for f in folders_to_delete:
            try:
                shutil.rmtree(f)
                print(f"Deleted: {f}")
            except Exception as e:
                print(f"Error deleting {f}: {e}")     
    