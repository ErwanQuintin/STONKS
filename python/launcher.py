'''
Created on 28 févr. 2023

@author: michel
'''
import os
from constants import PATHTO
from session import Session
from werkzeug.datastructures import FileStorage
    
if __name__ == '__main__':
    obsid = "0804670301"
    filepath = os.path.join(PATHTO.master_sources, f"P{obsid}EPX000OBSMLI0000.FIT")
    Session.clean_up(1)
    
    my_file = FileStorage(
        stream=open(filepath, "rb"),
        name=f"P{obsid}EPX000OBSMLI0000.FIT",
        content_type="octet/stream",
        )
    session = Session(my_file)    
    print(session.process_observation())
    session.get_tarball()

