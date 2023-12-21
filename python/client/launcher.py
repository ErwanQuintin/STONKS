'''
Created on 6 déc. 2023

@author: michel
'''
from client.list_processor import ListProcessor
if __name__ == '__main__':
    processor = ListProcessor()
    processor.mprocess_files(5)
