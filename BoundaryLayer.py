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


def slip_operator(index, Setax, Setay):
    Uzero = np.array([[Setay, -1*Setax], [Setax, Setay]])
    Uone = np.array([[Setay, -1*Setax], [-1*Setax, -1*Setay]])
    O = np.invert(Uone)*Uzero
    vecUzero = O*np.array([[state[index]], []])


def nonslip_operator():
    pass


def adiabatic():
    pass


def non_adiabatic():
    pass


def wall_operator(wall, energy, index, od, ob, ot):
    wall_cond = {"slip": nonslip_operator, "nonslip": slip_operator}
    energy = {"adiabatic": adiabatic, "nonadiabatic": non_adiabatic}

    # Instead of assuming geometry we take in centers with these conditions
    # Allow us to reuse this operator for any geometry
    # Backing out cell values at points
    
    yeta = totalh[od, 2] - totalh[ol, 2]
    xeta = totalh[od, 1] - totalh[ol, 1]
    yxi = totalh[od, 2] - totalh[ot, 2]
    xxi = totalh[od, 1] - totalh[ot, 1]
    
    Setax = -1*yxi
    Setay = xxi

    momentum_operator = wall_cond[wall](index, Setax, Setay)


# Initializing State Vector
Qi = state_vector(101325, 300, 2.000)  # Pa, K, Mach Numeber
size = np.shape(totalh)
state = np.zeros([size[0], size[1], len(Qi)])
state[3:-3:2, 1:-1:2] = Qi
for j in range(len(state[1::2, 1])):
    # Iterate over rows fixing the first column
    print(state[2*j+1, 1])
    state[:, 0][j] = wall_operator("slip", "adiabatic", [2*j + 1, 1])
