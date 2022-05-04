# -*- coding: utf-8 -*-
"""
Created on Mon May  2 22:21:39 2022

@author: Dhruva
"""
import numpy as np
from GridReader import reader, metrics, halo_augmenter



def read_in(file):
    data, grid, haloXL, haloXR, haloYT, haloYB = reader(file)
    gridh = halo_augmenter(grid, haloXL, haloXR, haloYB, haloYT)
    total = metrics(grid)
    totalh = metrics(gridh)
    
    return totalh

def state_vector(P, T, M, Rbar=287, c=347.2, gamma=1.4):
    # Density
    rho = Rbar * T/P
    q_one = rho

    # Velocity U
    u = 2*c
    q_two = rho*u

    # Velocity V
    v = 0
    q_three = rho*v

    # Energy
    density_energy = P/(gamma - 1) + 0.5*rho*(u**2 + v**2)
    q_four = density_energy

    Q = np.array([q_one, q_two, q_three, q_four])
    Qv = np.array([P, u, v, T])
    return Q, Qv

def density(index, P, T, state):
    #Assume normal pressure gradient at wall is 0, implying gradient = 0
    #Techinally rho = Po(Pin)/R*To(Tin)
    return state[index[0], index[1], 0]

def slip_operator(index, Setax, Setay, state):
    Uzero = np.array([[Setay, -1*Setax], [Setax, Setay]])
    Uone = np.array([[Setay, -1*Setax], [-1*Setax, -1*Setay]])
    Uizero = np.linalg.inv(Uzero)
    Op = np.dot(Uizero, Uone)
    vzero = np.array([state[index[0], index[1], 1], state[index[0], index[1], 2]])
    vecUzero = -1*np.dot(Op, vzero)
    
    return vecUzero


def nonslip_operator(x, y, z):
    pass


def adiabatic(index, P, state, gamma=1.4):
    total_energy = P/(gamma - 1) + 0.5*state[index[0], index[1], 0]*(state[index[0], index[1], 1]**2 + state[index[0], index[1], 2]**2)
    return total_energy


def non_adiabatic():
    pass


def wall_operator(wall, energy, index, totalh, state):
    wall_cond = {"nonslip": nonslip_operator, "slip": slip_operator}
    energy = {"adiabatic": adiabatic, "nonadiabatic": non_adiabatic}

    # Instead of assuming geometry we take in centers with these conditions
    # Allow us to reuse this operator for any geometry
    # Backing out cell values at points
    tl = [index[0] + 1, index[1] - 1]
    bl = [index[0] - 1, index[1] - 1]
    tr = [index[0] + 1, index[1] + 1]
    br = [index[0] - 1, index[1] + 1]
    yeta = totalh[tr[0], tr[1],1] - totalh[tl[0], tl[1], 1] 
    xeta = totalh[tr[0], tr[1], 0] - totalh[tl[0], tl[1], 0]
    yxi = totalh[tr[0], tr[1], 1] - totalh[br[0], br[1], 1]
    xxi = totalh[tr[0], tr[1], 0] - totalh[br[0], br[1], 0]
    
    Setax = -1*yxi
    Setay = xxi

    halo_momentum = wall_cond[wall](index, Setax, Setay, state)
    
    return halo_momentum
    

def initial_values(P, T, M, totalh):
    # Initializing State Vector
    Qi, Qv = state_vector(P, T, M)  # Pa, K, Mach Numeber
    size = np.shape(totalh)
    state = np.zeros([size[0], size[1], len(Qi)])
    state[3:-3:2, 1:-1:2] = Qi
    
    
    for j in range(len(state[3, 1::2])):
        # Iterate over rows fixing the first column
        index = [3, 2*j+1]
        rho = density(index, Qv[0], Qv[3], state)
        momentum = wall_operator("slip", "adiabatic", index, totalh, state)
        
        
        state[1, 2*j+1, 0] = rho
        state[1, 2*j + 1, 1] = momentum[0]
        state[1, 2*j + 1, 2] = momentum[1]
        
        energy = adiabatic(index, Qv[0], state)
        state[1, 2*j + 1, 3] = energy
        
        
        index = [-4, 2*j + 1]
        rho = density(index, Qv[0], Qv[3], state)
        momentum = wall_operator("slip", "adiabtic", index, totalh, state)
        energy = adiabatic(index, Qv[0], state)    
        
        state[-2, 2*j+1, 0] = rho
        state[-2, 2*j + 1, 1] = momentum[0]
        state[-2, 2*j + 1, 2] = momentum[1]
        
        energy = adiabatic(index, Qv[0], state)
        state[-2, 2*j + 1, 3] = energy
        
    return state
    