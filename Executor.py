# -*- coding: utf-8 -*-
"""
Created on Tue May  3 20:17:30 2022

@author: Dhruva
"""
import numpy as np
from BoundaryLayer import read_in, initial_values
from Solver import state_interpolation, Roe_Jacobian, Convective_Operator

M = 2 # Mach Number Axial
c = 347.2 # Sound Speed (m/s)
P = 101325 # Pressure (Pa)
T = 300 # Kelvin

file = "Coarse.dat"
totalh = read_in(file)
state = initial_values(P, T, M, totalh)

def time_step(totalh, state, c=347.2):
    size = np.shape(totalh)
    CFL_grid = np.zeros([size[0], size[1], 2])
    for i in range(3,size[0] - 3, 2):
        for j in range(3, size[1] - 3, 2):
            
            rho_xi = abs(state[i, j, 1]) + 0.5*c*(totalh[i+1, j, 1] + totalh[i-1, j, 1])
            rho_eta = abs(state[i, j, 2]) + 0.5*c*(totalh[i, j-1, 2] + totalh[i, j+1, 2])
            
            CFL_grid[i, j] = np.array([rho_xi, rho_eta])
    
    CFL = [np.min(CFL_grid[3:-3:2, 3:-3:2, 0]), np.min(CFL_grid[3:-3:2, 3:-3:2, 1])]
    data = CFL_grid
            
    return CFL, data
            

def iterator(state, totalh, max_iterations):
    delta_tau, data = time_step(totalh, state)
    p = 0.05
    delta_t = min(delta_tau) *(1 - p)
    
    size = np.shape(state)
    var = np.zeros([size[0], size[1], 2, 4])
    
    state = np.random.rand(size[0], size[1], 4)
    inp = state
    
    for i in range(3,size[0] - 3, 2):
        for j in range(3, size[1] - 3, 2):
            # print(state[i, j])
            # print(i, j)
            QB, QT = state_interpolation([i+1, j], 0, 1, "eta", state)
            var[i+1, j, 0] =  QB
            var[i+1, j, 1] = QT
                       
            QL, QR = state_interpolation([i, j+1], 0, 1, "xi", state)
            var[i, j+1, 0] = QL
            var[i, j+1, 1] = QR
            
            
        
    i=0
    Qprev = state
    # while i<max_iterations:
    #     Qnew = Qprev - delta_t*(2)
    
    return var, inp

x, inp =iterator(state, totalh, 500)
xvis = x[:, :, 1, :]