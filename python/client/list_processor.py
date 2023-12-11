'''
Created on 6 déc. 2023

@author: michel
'''
import os, sys
from urllib.parse import unquote
import requests
import tempfile
import multiprocessing as mp
import gzip

STONKS_URL = "http://serendib.astro.unistra.fr:5555/stonks/process-obs"
DATA_DIR = "/rawdata/4XMMdr13/ready_to_load/"

INPUT_DIR = os.path.join(DATA_DIR, "EpicSumSrcList")
OUTPUT_DIR = os.path.join("/data/STONKS", "EpicVarAlert")


class ListProcessor:
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        print(f"build the input file list")

        input_list = os.listdir(INPUT_DIR)
        self._input_files = []
        for _input in input_list:
            if 'OBSMLI' in _input:
                self._input_files.append(
                    os.path.realpath(os.path.join(INPUT_DIR, _input)))
        print(f"{len(self._input_files)} files to process")
    
    def _get_response_filepath(self, response):
        content_disposition = response.headers['Content-Disposition']

        # Extract the filename from the Content-Disposition header
        filename_index = content_disposition.find('filename=')
        if filename_index != -1:
            filename = content_disposition[filename_index + len('filename='):]
            filename = unquote(filename)  # Decode URL-encoded characters
            return os.path.join(OUTPUT_DIR, filename)
        else:
            return os.path.join(OUTPUT_DIR, "MissingFilename.tar")
        
    def _get_uncompress_filepath(self, file_path):
        if not file_path.endswith(".gz"):
            return file_path, False
        
        print(f"uncompress {file_path}...")
        temp_dir = tempfile.gettempdir() 
        out_file = os.path.join(temp_dir, os.path.basename(file_path)[:-3])

        with gzip.open(file_path, 'rb') as gzip_file:
            # Read the contents of the gzip file
            file_content = gzip_file.read()

        # Write the contents to a temporary file
        with open(out_file, 'wb') as temp_file:
            temp_file.write(file_content)
        print(f" ... in {out_file}")

        return out_file, True
    
    def _check_if_processed(self, obs_id):
        
        output_list = os.listdir(OUTPUT_DIR)
        for output in output_list:
            if obs_id in output:
                print(f"{obs_id} already processed")
                return True
        return False
    
    def _process_file(self, file_path):
        try:
            # Open the file to be sent
            print(f"process {file_path}")
            file_path, compressed = self._get_uncompress_filepath(file_path)
            obs_id = os.path.basename(file_path)[1:11]
            if self._check_if_processed(os.path.basename(file_path)[1:11]):
                return
            with open(file_path, 'rb') as file:
                # Prepare the file and any additional data if needed
                files = {'file': (file_path, file)}
                # Make the POST request to the service
                response = requests.post(STONKS_URL, files=files)
            
                print(f"Get the response from {STONKS_URL} ")
    
                # Check if the request was successful (status code 200)
                if response.status_code == 200:
                    # Save the response content to a file
                    output_file_path = self._get_response_filepath(response)
                    print(f"to be saved in {output_file_path}")
                    with open(output_file_path, 'wb') as output_file:
                        output_file.write(response.content)
                    print(f"File successfully sent and response saved to {output_file_path}")
                else:
                    print(f"Error: Received status code {response.status_code} from the service {response.content}")
                    with open(os.path.join(OUTPUT_DIR, f"{obs_id}.empty"), 'wb') as output_file:
                        output_file.write(b"empty file\n")

                    if response.status_code != 404:
                        sys.exit(1)

            if compressed:
                print(f"remove {file_path}")
                os.remove(file_path)
        except Exception as e:
            print(f"Error: {e}")

    def process_files(self):
        
        for file_path in self._input_files:
            self._process_file(file_path)
            break

    def mprocess_files(self, nb_process):
        
        with mp.Pool(nb_process) as pool:
            pool.map(self._process_file, self._input_files)

        