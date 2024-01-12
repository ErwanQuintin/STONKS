'''
Created on 1 mars 2023

@author: michel
'''
import os

class PATHTO:
    basedir = os.path.join(os.path.dirname(__file__), "..")
    catalogs = os.path.join(basedir, "Data/Catalogs/FullData/")
    master_sources =  os.path.join(basedir, "Data/MasterSource/")
    stilts_cmd = os.path.join(basedir, "stilts/stilts")
    sessions = os.path.join(basedir, 'sessions')
    precomputed_obsids = os.path.join(master_sources, 'PreComputedObsidMatches')
print(f"working from {PATHTO.sessions}")
class CACHE:
    precomputed_depth = 50
    number_of_cached_sessions = 1000
    cleanup_delay = 5000 # in sec

#if os.path.exists(PATHTO.tempo) is False:
#    print(f"create tempo dir {PATHTO.tempo}")
#    os.mkdir(PATHTO.tempo)
