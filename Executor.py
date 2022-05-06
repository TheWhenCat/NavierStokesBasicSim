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
    
    A_average = np.zeros([size[0], size[1], 4, 4])
    E = np.zeros([size[0], size[1], 4])
    print(np.shape(E))
    F = np.zeros([size[0], size[1], 4])
    m= 3
    xmax = size[0] - 3
    for i in range(m,xmax, 2):
        k = 3
        c=0
        maximum = size[1] - 3
        for j in range(k, maximum, 2):
            if all(inp[i, j] == 0):
                print(inp[i, j])
                
            QB,QT = state_interpolation([i, j], 0, 1, "eta", state)
            # logging.debug(f'{QB}')
            # logging.debug(f'{QT}')
            var[i-1, j, 0] = QB
            var[i+1, j, 1] = QT
                       
            QL, QR = state_interpolation([i, j+1], 0, 1, "xi", state)
            var[i, j-1, 0] = QL
            var[i, j+1, 1] = QR
            
            #normalize    
            neta = (totalh[i-1, j-1, 0:2] - totalh[i-1, j+1, 0:2])
            neta = neta / np.sqrt(neta.dot(neta))
            nxi = (totalh[i+1, j-1, 0:2] - totalh[i-1, j+1, 0:2])
            print(nxi)
            nxi = nxi / np.sqrt(nxi.dot(nxi))
            try:
                A_average[i+1, j] = Roe_Jacobian(neta, QB, QT)
                A_average[i, j+1] = Roe_Jacobian(nxi, QL, QR)
            except:
                print([i, j])
            print(Convective_Operator(nxi, state[i+1, j]))
            E[i+1, j] = np.array(Convective_Operator(nxi, state[i+1, j]))
            F[i, j+1] = np.array(Convective_Operator(neta, state[i, j+1]))
            
            if c<6:
                c+=1
            else:
                break
        #     k = k + 2
        #     if k+1 == maximum:
        #         break
        # m = m + 2
        # if m + 1 == xmax:
        #     break
        
    i=0
    Qprev = state
    while i<0: #max_iterations:
        Epos2 = 0 #np.matmul(E[:: ,:], A_average[:, 1:])
        Spos2 = 0
        Eneg2 = 0
        Sneg2 = 0
        
        Fpos2 = 0
        Fneg2 = 0
        
        
        Qnew = Qprev - delta_t*( (Epos2 * Spos2 - Eneg2 * Sneg2 ) + (Fpos2 * Spos2 - Fneg2 * Sneg2))
    
    return var, E, F

