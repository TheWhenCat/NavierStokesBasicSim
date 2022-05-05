# -*- coding: utf-8 -*-
"""
Created on Tue May  3 20:28:07 2022

@author: Dhruva
"""
import numpy as np

def Roe_Jacobian(eta, stateL, stateR, gamma = 1.4):
    #Want to operate using a Roe Averages Scheme of the Values
    #hat{A}
    
    etax = eta[0]
    etay = eta[1]
    
    rhoL = stateL[0]
    uL = stateL[1]
    vL = stateL[2]
    hL = stateL[3]
    
    rhoR = stateR[0]
    uR = stateR[1]
    vR = stateR[2]
    hR = stateR[3]
    
    rho = np.sqrt(rhoL * rhoR)
    u = (np.sqrt(rhoR) * uR + np.sqrt(rhoL) * uL) *1/(np.sqrt(rhoR) + np.sqrt(rhoL))
    v = (np.sqrt(rhoR) * vR + np.sqrt(rhoL) * vL) *1/(np.sqrt(rhoR) + np.sqrt(rhoL))
    h = (np.sqrt(rhoR) * hR + np.sqrt(rhoL) * hL) * 1/(np.sqrt(rhoR) + np.sqrt(rhoL))
    c = 
    
    
def state_interpolation():
    pass


