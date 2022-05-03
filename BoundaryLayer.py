# -*- coding: utf-8 -*-
"""
Created on Mon May  2 22:21:39 2022

@author: Dhruva
"""
import numpy as np
from GridReader import reader, metrics, halo_augmenter

file = "Coarse.dat"
data, grid, haloXL, haloXR, haloYT, haloYB = reader(file)
gridh = halo_augmenter(grid, haloXL, haloXR, haloYB, haloYT)
total = metrics(grid)
totalh = metrics(gridh)

def state_vector(P, T, M, Rbar=287, c=347.2, gamma=1.4):
    #Density
    rho = Rbar * T/P
    q_one = rho
    
    #Velocity U
    u = 2*c
    q_two = rho*u
    
    #Velocity V
    v = 0
    q_three = rho*v
    
    #Energy
    density_energy = P/(gamma - 1) + 0.5*rho*(u**2 + v**2)
    q_four = density_energy
    
    Q = np.array([q_one, q_two, q_three, q_four])
    return Q

def wall_operator(wall, energy, center, grid, total):
    #Instead of assuming geometry we take in centers with these conditions
    #Allow us to reuse this operator for any geometry
    yeta = grid
    xeta =
    yxi = 
    xxi = 
    

#Initializing State Vector
Qi = state_vector(101325, 300, 2.000)
size = np.shape(totalh)
state = np.empty([size[0], size[1], len(Qi)])
state[:, :] = Qi
for i in grid[0]


