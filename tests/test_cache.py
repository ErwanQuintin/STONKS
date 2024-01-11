'''
Created on 18 Dec 2023

@author: laurentmichel
'''
import os
import unittest
import time
import glob
import shutil
from cache import Cache

data_dir = os.path.join(os.path.dirname(__file__), "data", "cache")

class Test(unittest.TestCase):


    def testPrecomputed(self):
        obsid = 1234567890
        shutil.rmtree(data_dir)
        os.mkdir(data_dir)

        for _ in range(40): 
            print(f"create files for upperlimit {obsid} ")            
            with open(f"{data_dir}/{obsid}.fits", mode='a'): pass
            with open(f"{data_dir}/UpperLimit{obsid}.fits", mode='a'): pass
            obsid += 1
            time.sleep(1)

        Cache._clean_precomputed(data_dir, 12, day=1)
        self.assertEqual(len(glob.glob(data_dir + '/*.fits', recursive=False)), 18)      
        Cache._clean_precomputed(data_dir, 12, day=1)
        self.assertEqual(len(glob.glob(data_dir + '/*.fits', recursive=False)), 18) 
        contents = os.listdir(data_dir)
        print(glob.glob(data_dir + '/*.fits', recursive=False))
        # Iterate through each item in the directory
        shutil.rmtree(data_dir)
        os.mkdir(data_dir)
             
    def testSessions(self):
        shutil.rmtree(data_dir)
        os.mkdir(data_dir)

        obsid = 1234567890
        for _ in range(40): 
            print(f"create session for obs {obsid} ")            
            os.mkdir(f"{data_dir}/Session_{obsid}") 

            with open(f"{data_dir}/Session_{obsid}/output_{obsid}.fits", mode='a'): pass
            obsid += 1
            time.sleep(1)
            
        Cache._clean_sessions(data_dir, 12, day=1)
        self.assertEqual(len(glob.glob(data_dir + '/*', recursive=False)), 9)      
        Cache._clean_sessions(data_dir, 12, day=1)
        self.assertEqual(len(glob.glob(data_dir + '/*', recursive=False)),9)      
        contents = os.listdir(data_dir)
        print(glob.glob(data_dir + '/Session*', recursive=False))
        # Iterate through each item in the directory
        shutil.rmtree(data_dir)
        os.mkdir(data_dir)

        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
