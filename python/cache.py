'''
Created on 1 mars 2023

@author: michel
'''
import os
import re
import shutil
import glob
import time
from threading import Thread
from constants import PATHTO, CACHE

class Cache(Thread):
    running = True
    
    @staticmethod
    def _clean_precomputed(cache_dir, max_length, day=86400):
        print(f"cleaning {cache_dir}: size limited to {max_length} obs")
        files = glob.glob(cache_dir + '/[0-9]*', recursive=False)      
        if len(files)  < max_length:
            print(f"{len(files)}  precomputed corrs in cache: not need to clean up") 
            return
        # get the current time 
        current_time = time.time()

        # loop over all the files 
        for ftr in files:
            file_time = os.stat(ftr).st_mtime
            dt = current_time - 10*day
            if(file_time < dt):
                obs = re.search("^.*([0-9]{10}).*$", ftr).group(1)
                files = glob.glob(cache_dir + f'/*{obs}*', recursive=False)    
                for _ftr in files:
                    os.remove(_ftr)  
                    print(f" Delete : {_ftr} dt={(current_time - file_time)}")
        return

    @staticmethod
    def _clean_sessions(session_dir, max_length, day=86400):
        
        print(f"cleaning {session_dir}: size limited to {max_length} obs")
        files = glob.glob(session_dir + '/*', recursive=False)      
        # get the current time 
        current_time = time.time() 
        if len(files)  < max_length:
            print(f" {len(files)}  sessions in cache: not need to clean up") 
            return
        # loop over all the files 
        for session in files: 
            file_time = os.stat(session).st_mtime 
            dt = current_time - 10*day
            if(file_time < dt): 
                print(f" Delete : {session}  dt={(current_time - file_time)} ") 
                if os.path.isdir(session):
                    shutil.rmtree(session)
                else:
                    os.remove(session)

    def run(self):
        print("Cache cleaner thread started")
        while Cache.running:
            print("Cache cleanup")
            Cache._clean_precomputed(PATHTO.precomputed_obsids, CACHE.precomputed_depth)
            Cache._clean_sessions(PATHTO.sessions, CACHE.number_of_cached_sessions)
            time.sleep(CACHE.cleanup_delay)
              
    
