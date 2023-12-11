'''
Created on 28 févr. 2023

@author: michel
'''
import os
from constants import PATHTO
from session import Session
    
if __name__ == '__main__':
    obsid = "0804670301"
    filepath = os.path.join(PATHTO.master_sources, f"P{obsid}EPX000OBSMLI0000.FIT")
    Session.clean_up(1)
    session = Session(filepath)
    session.process_observation(None)
    session.get_tarball()
