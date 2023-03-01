'''
Created on 28 févr. 2023

@author: michel
'''
import os
from rest_api.logic import process_one_observation
from constants import PATHTO
from session_utils import SessionUtils

    
if __name__ == '__main__':
    obsid = "0804670301"
    SessionUtils.remove_session_directory(obsid)
    process_one_observation(os.path.join(PATHTO.master_sources, f"P{obsid}EPX000OBSMLI0000.FIT"))
