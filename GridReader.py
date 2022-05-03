# -*- coding: utf-8 -*-
"""
Created on Sun Mar 20 17:39:20 2022

@author: Dhruva
"""

# import pandas as pd
import numpy as np
import itertools
# from Hashtable import Hashtable
import hashlib


def is_float(string):
    """ True if given string is float else False"""
    try:
        return float(string)
    except ValueError:
        return False
    
def combine(x,y):
    val = (2**x)*(2*y + 1) - 1
    return val
    
def Convert(tup, di):
    for a, b in tup:
        di.setdefault(a, []).append(b)
    return di

def reader(file):
    data = []
    #Reading the file
    file = "Coarse.dat"
    with open(file, 'r') as f:
        d = f.readlines()[1:]
        for i in d:
            k = i.rstrip().split(",")
            data.append([float(i) for i in k]) 
    
    #Find the discrete values
    data = np.array(data, dtype='O')
    combined = [(hashlib.sha256(str((i[0], i[1])).encode('utf-8')).hexdigest(), (i[0], i[1])) for i in data]
    
    hash_proper = {}
    hash_proper = Convert(combined, hash_proper)
    X = np.unique(data[:, 0])
    Y = np.unique(data[:, 1])
    
    grid = []
    for i in range(len(X)):
        row = []
        for j in range(len(Y)):        
            val = hashlib.sha256(str((X[i], Y[j])).encode('utf-8')).hexdigest()
            if val in hash_proper.keys():
                    row.append(np.array([X[i], Y[j]]))
        grid.append(np.array(row))
        
    grid = np.array(grid)
    
    haloXL = np.array([grid[0, :, 0] - (grid[1, :, 0] - grid[0, :, 0]), grid[0, :, 1]])
    haloXR = np.array([grid[-1, :, 0] - (grid[-2, :, 0] - grid[-1, :, 0]), grid[-1, :, 1]])
    haloYT = np.array([grid[:, -1, 0], grid[:, -1, 1] + (grid[:, -1, 1] - grid[:, -2, 1])])
    haloYB = np.array([grid[:, 0, 0], grid[:, 0, 1] + (grid[:, -2, 1] - grid[:, -1, 1])])
    
    haloXR = np.insert(haloXR, len(haloXR[0]),  [grid[-1, -1, 0] - (grid[-2, -1, 0] - grid[-1, -1, 0]),  grid[-1, -1, 1] - (grid[-1, -2, 1] - grid[-1, -1, 1])], axis=1)
    haloXR = np.insert(haloXR, 0, [grid[-1, 0, 0] - (grid[-2, 0, 0] - grid[-1, 0, 0]),  grid[-1, 0, 1] - (grid[-1, 1, 1] - grid[-1, 0, 1])], axis = -1)
    
    haloXL = np.insert(haloXL, len(haloXL[0]),  [grid[0, -1, 0] - (grid[1, -1, 0] - grid[0, -1, 0]),  grid[0, -1, 1] - (grid[1, -2, 1] - grid[0, -1, 1])], axis=1)
    haloXL = np.insert(haloXL, 0, [grid[0, 0, 0] - (grid[1, 0, 0] - grid[0, 0, 0]),  grid[0, 0, 1] - (grid[0, 1, 1] - grid[0, 0, 1])], axis = -1)
    
    return data, grid, haloXL, haloXR, haloYT, haloYB

def metrics(grid):
    #Define an array to store values with its associated coordinate
    total = np.zeros([len(grid)*2- 1, len(grid[0])*2 -1, 3])
    #Inserting coordinates at even, even pairs
    size = np.shape(total)
    for i in range(0:size[0]:2):
        for j in range(0:size[1]:2):
            total[i, j, 1:2] = grid[i/2, j/2]
            
    #Volume Calculation
    for i in range(len(grid) - 1):
        l = len(grid[i])
        for j in range(l-1):
            cen =  [0.25*(grid[i][j][0] + grid[i+1][j][0] + grid[i][j+1][0] + grid[i+1][j+1][0]), 0.25*(grid[i][j][1] + grid[i+1][j][1] + grid[i][j+1][1] + grid[i+1][j+1][1])]
            vol = .5*((grid[i+1][j+1][0] - grid[i][j][0])*(grid[i][j+1][1] - grid[i+1][j][1]) -  (grid[i][j+1][0] - grid[i+1][j][0])*(grid[i+1][j+1][1] - grid[i][j][1]))
            cen.append(vol)
            #Any possible two odds will store the cell volumes
            total[2*i + 1, 2*j + 1] = cen
            
    #Flux 
    for i in range(len(grid) - 1):
        for j in range(len(grid[i]) - 1):
            #Si fluxes
            fc = []
            f = []
            fc = [0.5*(grid[i+1][j][0] + grid[i][j][0]),  0.5*(grid[i][j][1] + grid[i+1][j][1])]
            f = np.sqrt((grid[i+1][j][0] - grid[i][j][0])**2 + (grid[i+1][j][1] - grid[i][j][1])**2)
            fc.append(f)
            #Vertical Fluxes in Eta direction
            #Stored in odd rows and even columns
            total[2*i + 1, 2*j ] = fc
            
            #Eta Fluxes
            fc = []
            f = []
            fc = [0.5*(grid[i][j][0] + grid[i][j+1][0]), 0.5*(grid[i][j][1]+ grid[i][j+1][1])]
            f = np.sqrt((grid[i][j][0] - grid[i][j+1][0])**2 + (grid[i][j][1]- grid[i][j+1][1])**2)
            fc.append(f)
            #Horizontal Fluxes in Xi direction
            #Stored in even rows and odd columns
            total[2*i , 2*j + 1] = fc
            
    for j in range(len(grid[-1, :]) - 1):
        fc = [grid[-1, j, 0], 0.5*(grid[-1,j,1] + grid[-1, j+1, 1])]
        f = -1*(grid[-1,j,1] - grid[-1, j+1, 1])
        fc.append(f)
        total[-1, 2*j+1] = fc
            
    for i in range(len(grid[:, -1]) - 1):
        fc = [0.5*(grid[i, -1, 0] + grid[i+1, -1, 0]), grid[i, -1, 1]]
        f = -1*(grid[i, -1, 0] - grid[i+1, -1, 0])
        fc.append(f)
        total[2*i + 1, -1] = fc
        
    return total
    
def halo_augmenter(grid, XL, XR, YB, YT):
    size = np.shape(grid)
    gridh = np.zeros([size[0]+2, size[1]+2, size[2]])
    gridh[1:-1, 0] = np.transpose(YB)
    gridh[1:-1, -1] = np.transpose(YT)
    gridh[0] = np.transpose(XL)
    gridh[-1] = np.transpose(XR)
    
    gridh[1:-1, 1:-1] = grid
    
    return gridh
    
    
    





    
