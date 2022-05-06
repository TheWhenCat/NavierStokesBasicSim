# -*- coding: utf-8 -*-
"""
Created on Tue May  3 20:17:30 2022

@author: Dhruva
"""
import time
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
    # print(np.shape(E))
    F = np.zeros([size[0], size[1], 4])
    m= 3
    xmax = size[0] - 4
    for i in range(m,xmax, 2):
        k = 3
        maximum = size[1] - 4
        for j in range(k, maximum, 2):
            if all(inp[i, j] == 0):
                print(inp[i, j])
                
            QB,QT = state_interpolation([i, j], 1, -1, "eta", state)
            # print("QB: ", QB)
            # print("QT: ", QT)
            # print("Difference: ", np.subtract(QB, QT))
            var[i-1, j, 0] = QB
            var[i+1, j, 1] = QT
                       
            QL, QR = state_interpolation([i, j], 1, -1, "xi", state)
            # print("QL: ", QL)
            # print("QR: ", QR)
            # print("Difference: ", np.subtract(QR, QL))
            var[i, j-1, 0] = QL
            var[i, j+1, 1] = QR
            
            #normalize    
            neta = (totalh[i-1, j-1, 0:2] - totalh[i-1, j+1, 0:2])
            neta = neta / np.sqrt(neta.dot(neta))
            nxi = (totalh[i+1, j-1, 0:2] - totalh[i-1, j+1, 0:2])
            
            nxi = nxi / np.sqrt(nxi.dot(nxi))
            A_average[i, j+1] = Roe_Jacobian(neta, QB, QT)
            A_average[i+1, j] = Roe_Jacobian(nxi, QL, QR)
            #print("Average: ", A_average[i+1, j])
            
            one = 0.5*(np.array(Convective_Operator(nxi, QL)) + np.array(Convective_Operator(nxi, QR)))
            two = -0.5*A_average[i+1, j].dot(np.subtract(QL, QR))
            # print("One: ", one)
            # print("Two: ", two)
            E[i+1, j] =   (one*state[i+1, j, 2] + two) * totalh[i+1, j, 2]
            # print("E: ", E[i+1, j])
            
            one = 0.5*(np.array(Convective_Operator(neta, QB)) + np.array(Convective_Operator(neta, QT)))
            two = -0.5*A_average[i, j+1].dot(np.subtract(QB, QT))
            F[i, j+1] = (one*state[i, j+1, 2] + two) * totalh[i, j+1, 2]
            # print("F: ", F[i, j+1])
            
        return var, E, F, delta_t
    
def total_flux(totalh, E, F, delta_t, Qprev, maxi=5):
    Etot = np.zeros(np.shape(Qprev)) 
    Ftot = np.zeros(np.shape(Qprev))
    i = 0
    t0 = time.time()
    while i<maxi: #max_iterations:
        # print("Everything is not ok :(")
        Etot[3::2, 0::2] = E[3::2, 0::2]- E[1:-3:2, 0::2]
        
        Ftot[3::2, 0::2] = F[3::2, 0::2] - F[1:-3:2, 0::2]
        
        Qnew = Qprev - delta_t*( Etot + Ftot)
        #res = np.sqrt(Qnew - Qprev).dot(Qnew - Qprev))
        # print("JK we made it")
        i += 1
        
        Qres = Qnew - Qprev
        sum = np.zeros(np.shape(Qres))
        for i in range(len(Qres)):
            for j in range(len(Qres[i])):
                sum += np.square(Qres[i, j])
        res = np.sqrt(sum)
    print("Time: ", time.time() - t0)
    # print("Residual: ", res)
    return Qnew #, res
    

