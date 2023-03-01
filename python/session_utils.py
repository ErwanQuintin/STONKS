'''
Created on 1 mars 2023

@author: michel
'''
import os
import tarfile

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
 
    @staticmethod
    def get_session_path(session_name):      
        return os.path.realpath(os.path.join(PATHTO.alert_light_curves,
                            session_name))
     
     
    @staticmethod
    def get_tarball(session_name):      
        output = os.path.join(PATHTO.tempo, f"{session_name}.tgz")
        sourcedir = SessionUtils.get_session_path(session_name)
        with tarfile.open(output, "w:gz") as tar:
            tar.add(sourcedir, arcname=os.path.basename(sourcedir))
        return output
   