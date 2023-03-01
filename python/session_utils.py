'''
Created on 1 mars 2023

@author: michel
'''
import os
import shutil
from constants import PATHTO

class SessionUtils(object):
    '''
    Some utilities helping to manage the session data
    '''
    @staticmethod
    def remove_session_directory(session_name):

        for root, dirs, files in os.walk(PATHTO.alert_light_curves):
            for f in files:
                os.unlink(os.path.join(root, f))
            for d in dirs:
                if d == session_name:
                    shutil.rmtree(os.path.join(root, d))
 
        