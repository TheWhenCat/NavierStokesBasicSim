# -*- coding: utf-8 -*-
"""
Created on Tue May  3 20:17:30 2022

@author: Dhruva
"""
import logging, sys
logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)

import numpy as np
from Solver import state_interpolation, Roe_Jacobian, Convective_Operator


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
    
    #state = np.random.rand(size[0], size[1], 4)
    inp = state
    A_average = np.zeros([size[0], size[1], 4])
    E = np.zeros([size[0], size[1], 4])
    F = np.zeros([size[0], size[1], 4])
    for i in range(3,size[0] - 3, 2):
        for j in range(3, size[1] - 3, 2):
            
            QB,QT = state_interpolation([i, j], 0, 1, "eta", state)
            logging.debug(f'{QB}')
            logging.debug(f'{QT}')
            var[i-1, j, 0] = QB
            var[i+1, j, 1] = QT
                       
            QL, QR = state_interpolation([i, j+1], 0, 1, "xi", state)
            var[i, j-1, 0] = QL
            var[i, j+1, 1] = QR
            
            #normalize    
            neta = (totalh[i-1, j-1, 0:2] - totalh[i-1, j+1, 0:2]) / np.sqrt(totalh[i-1, j-1, 0:2] - totalh[i-1, j+1, 0:2])
            nxi = (totalh[i+1, j-1, 0:2] - totalh[i-1, j+1, 0:2]) / np.sqrt(totalh[i+1, j-1, 0:2] - totalh[i-1, j+1, 0:2])
            
            A_average[i+1, j] = Roe_Jacobian(neta, QB, QT)
            A_average[i, j+1] = Roe_Jacobian(nxi, QL, QR)
            
            E[i+1, j] = Convective_Operator(nxi, state[i+1, j])
            F[i, j+1] = Convective_Operator(neta, state[i, j+1])
            
            
    i=0
    Qprev = state
    # while i<max_iterations:
    #     Qnew = Qprev - delta_t*(2)
    
    return var, inp

