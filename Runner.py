# -*- coding: utf-8 -*-
"""
Created on Thu May  5 06:18:20 2022

@author: Dhruva
"""
from Executor import iterator
from BoundaryLayer import read_in, initial_values

M = 2 # Mach Number Axial
c = 347.2 # Sound Speed (m/s)
P = 101325 # Pressure (Pa)
T = 300 # Kelvin

file = "Coarse.dat"
totalh = read_in(file)
state = initial_values(P, T, M, totalh)

x, inp =iterator(state, totalh, 500)
xvis = x[:, :, 0, 3]