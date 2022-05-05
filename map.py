# -*- coding: utf-8 -*-
"""
Created on Thu May  5 03:12:36 2022


"""

def N(data, i, j, k=0, skipby=1):
    if (k!=0):    
        return data[skipby*i-1][j][k]
    else:
        return data[skipby*i-1][j]

def S(data, i, j, k=0, skipby=1):
    if (k!=0):  
        return data[skipby*i+1][j][k]
    else:
        return data[skipby*i+1][j]

def E(data, i, j, k=0, skipby=1):
    if (k!=0):  
        return data[i][skipby*j+1][k]
    else:
        return data[i][skipby*j+1]

def W(data, i, j, k=0, skipby=1):
    if (k!=0):  
        return data[i][skipby*j-1][k]
    else:
        return data[i][skipby*j-1]

def NE(data, i, j, k=0, skipby=1):
    if (k!=0):  
        return data[i-1][skipby*j+1][k]
    else:
        return data[i-1][skipby*j+1]

def NW(data, i, j, k=0, skipby=1):
    if (k!=0):  
        return data[skipby*i-1][skipby*j-1][k]
    else:
        return data[skipby*i-1][skipby*j-1]
        
    
def SE(data, i, j, k=0, skipby=1):
    if (k!=0):  
        return data[skipby*i+1][skipby*j+1][k]
    else:
        return data[skipby*i+1][skipby*j+1]

def SW(data, i, j, k=0, skipby=1):
    if (k!=0):  
        return data[skipby*i+1][skipby*j-1][k]
    else:
        return data[skipby*i+1][skipby*j-1]

def cen(data, i, j, k=0, skipby=1):
    if (k!=0):  
        return data[i][j][k]
    else:
        return data[i][j]

