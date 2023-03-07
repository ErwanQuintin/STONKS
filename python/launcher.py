'''
Created on 28 févr. 2023

@author: michel
'''
import os
from rest_api.logic import process_one_observation
from constants import PATHTO
from session import Session
    
if __name__ == '__main__':
    obsid = "0804670301"
    Session.clean_up(1)
    session = Session()
    process_one_observation(session, os.path.join(PATHTO.master_sources, f"P{obsid}EPX000OBSMLI0000.FIT"),  None)
    session.get_tarball()
