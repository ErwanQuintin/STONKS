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
        for _ in range(40): 
            print(f"create files for obs {obsid} ")            
            with open(f"{data_dir}/{obsid}.fits", mode='a'): pass
            with open(f"{data_dir}/UpperLimit{obsid}.fits", mode='a'): pass
            obsid += 1
            time.sleep(1)

        Cache._clean_precomputed(data_dir, 12)
        self.assertEqual(len(glob.glob(data_dir + '/[0-9]*', recursive=False)), 12)      
        Cache._clean_precomputed(data_dir, 12)
        self.assertEqual(len(glob.glob(data_dir + '/[0-9]*', recursive=False)), 12) 
        contents = os.listdir(data_dir)

        # Iterate through each item in the directory
        for item in contents:
            item_path = os.path.join(data_dir, item)
            if os.path.isfile(item_path):
                os.remove(item_path)
             
    def testSessions(self):
        obsid = 1234567890
        for _ in range(40): 
            print(f"create files for obs {obsid} ")            
            with open(f"{data_dir}/{obsid}.fits", mode='a'): pass
            with open(f"{data_dir}/UpperLimit{obsid}.fits", mode='a'): pass
            obsid += 1
            time.sleep(1)
            
        Cache._clean_sessions(data_dir, 12)
        self.assertEqual(len(glob.glob(data_dir + '/*', recursive=False)), 12)      
        Cache._clean_sessions(data_dir, 12)
        self.assertEqual(len(glob.glob(data_dir + '/*', recursive=False)), 12)      
        contents = os.listdir(data_dir)
        # Iterate through each item in the directory
        for item in contents:
            item_path = os.path.join(data_dir, item)
            if os.path.isfile(item_path):
                os.remove(item_path)

        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()