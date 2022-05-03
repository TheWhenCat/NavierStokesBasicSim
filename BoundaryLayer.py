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
    return Q


def slip_operator():
    pass
    # Uzero =
    # Uone =


def nonslip_operator():
    pass


def adiabatic():
    pass


def non_adiabatic():
    pass


def wall_operator(wall, energy, index, orientation="None"):
    wall_cond = {"slip": nonslip_operator, "nonslip": slip_operator}
    energy = {"adiabatic": adiabatic, "nonadiabatic": non_adiabatic}
    directions = {"top": [[index[0] - 1, index[1]+1], [index[0]+1, index[1]+1]], 
                  "bottom": [[index[0] - 1, index[1]-1], [index[0]+1, index[1]-1]],
                  "left": [[index[0]-1, index[1]-1], [index[0]-1, index[1]+1]], 
                  "right": [[index[0]+1, index[1] - 1], [index[0]+1, index[1]+1]]
                  }
    # Instead of assuming geometry we take in centers with these conditions
    # Allow us to reuse this operator for any geometry
    # Backing out cell values at points
    ij = []
    tr = 
    yeta = 
    xeta = 0
    yxi = 0
    xxi = 0

    momentum_operator = wall_cond[wall]


# Initializing State Vector
Qi = state_vector(101325, 300, 2.000)  # Pa, K, Mach Numeber
size = np.shape(totalh)
state = np.zeros([size[0], size[1], len(Qi)])
state[1::2, 1::2] = Qi
for j in range(len(state[1::2, 1])):
    # Iterate over rows fixing the first column
    print(state[2*j+1, 1])
    state[:, 0][j] = wall_operator("slip", "adiabatic", [2*j + 1, 1])
