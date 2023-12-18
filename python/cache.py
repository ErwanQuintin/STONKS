'''
Created on 1 mars 2023

@author: michel
'''
import os
import re
import shutil
import glob
from threading import Thread
from constants import PATHTO, CACHE

class Cache(Thread):
    running = True
    
    @staticmethod
    def _clean_precomputed(cache_dir, max_length):
        print(f"cleaning {cache_dir}: size limited to {max_length} obs")
        files = glob.glob(cache_dir + '/[0-9]*', recursive=False)      

        files.sort(key=os.path.getctime)
        to_remove = files[max_length:]
        obs_to_remove = []
        for file in to_remove:
            obs_to_remove.append(re.search("^.*([0-9]{10}).*$", file).group(1))
        print(f"{len(obs_to_remove)} obs to remove")
        for obs in obs_to_remove:
            print(obs)
            files = glob.glob(cache_dir + f'/*{obs}*', recursive=False)    
            for ftr in files:
                os.remove(ftr)   
                
    @staticmethod
    def _clean_sessions(session_dir, max_length):
        print(f"cleaning {session_dir}: size limited to {max_length} obs")
        files = glob.glob(session_dir + '/*', recursive=False)      

        files.sort(key=os.path.getctime)
        session_to_remove = files[max_length:]
        print(f"{len(session_to_remove)} sessions to remove")
        for session in session_to_remove:
            if os.path.isdir(session):
                shutil.rmtree(session)  
            else:
                os.remove(session)   

    def run(self):
        while Cache.running:
            print("Cache cleanup")
            Cache._clean_precomputed(PATHTO.precomputed_obsids, CACHE.precomputed_depth)
            Cache._clean_sessions(PATHTO.sessions, CACHE.number_of_cached_sessions)
            self.sleep(CACHE.cleanup_delay)
              
    