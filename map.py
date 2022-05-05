# -*- coding: utf-8 -*-
"""
Created on Thu May  5 03:12:36 2022


"""

def N(data, i, j, k=0):
    return data[i-1][j][k]

def S(data, i, j, k=0):
    return data[i+1][j][k]

def E(data, i, j, k=0):
    return data[i][j+1][k]

def W(data, i, j, k=0):
    return data[i][j-1][k]

def NE(data, i, j, k=0):
    return data[i-1][j+1][k]

def NW(data, i, j, k=0):
    return data[i-1][j-1][k]

def SE(data, i, j, k=0):
    return data[i+1][j+1][k]

def SW(data, i, j, k=0):
    return data[i+1][j-1][k]

def cen(data, i, j, k=0):
    return data[i][j][k]

