# -*- coding: utf-8 -*-
"""
Created on Sun May  1 17:36:36 2022

@author: Dhruva
"""
import pprint

class Hashtable:
    def __init__(self, elements):
        self.bucket_size = len(elements) * 5
        self.buckets = [[] for i in range(self.bucket_size)]
        self._assign_buckets(elements)
        
    def _assign_buckets(self, elements):
        print("FUCKIUNG GET HERE")
        for key, value in elements:
            hashed_value = hash(key)
            index = hashed_value % self.bucket_size
            self.buckets[index].append((key, value))
            
    def get_value(self, input_key):
        hashed_value = hash(input_key)
        index = hashed_value % self.bucket_size
        bucket = self.buckets[index]
        for key,value in bucket:
            if key == input_key:
                return(value)
        return None
    
    def __str__(self):
        return pprint.pformat(self.buckets)
    

if __name__ == "__main__":
    capitals = [
        ('France', 'Paris'),
        ('United States', 'Washington D.C.'),
        ('Italy', 'Rome'),
        ('Canada', 'Ottawa')
        ]
    hashtable = Hashtable(capitals)
    print(hashtable)
    print(f"The capital is {hashtable.get_value('Italy')}")