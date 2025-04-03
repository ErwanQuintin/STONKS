'''
Created on 28 févr. 2023

@author: michel
'''
import os
from constants import PATHTO
from session import Session
from werkzeug.datastructures import FileStorage
    
if __name__ == '__main__':
    for varalert in (1,2,3):
        name = f'FAKE_0693760101_VARALERT{varalert}.gz'
        filepath = os.path.join(PATHTO.master_sources, name)
        #Session.clean_up(1)

        my_file = FileStorage(
            stream=open(filepath, "rb"),
            name=name,
            content_type="octet/stream",
            )
        session = Session(my_file)
        print(session.process_observation())
        session.get_tarball()

