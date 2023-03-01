'''
Created on 1 mars 2023

@author: michel
'''
from builtins import staticmethod

class DictUtils:
    """
    A few utilities facilitating the dictionary use here and there
    """
    @staticmethod
    def get_value_by_key(input_dict, key):
        """
        return ether the value attached to the key or "----"
        """
        if key in input_dict:
            return input_dict[key]
        return "----"
        
        