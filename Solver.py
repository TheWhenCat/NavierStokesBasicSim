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
    uL = stateL[1]/rhoL
    vL = stateL[2]/rhoL
    hL = (gamma-1)*rhoL*(stateL[3] - 0.5*(uL**2 + vL**2))
    
    rhoR = stateR[0]
    uR = stateR[1]/rhoR
    vR = stateR[2]/rhoR
    hR = (gamma-1)*rhoR*(stateR[3] - 0.5*(uR**2 + vR**2))
    
    rho = np.sqrt(rhoL * rhoR)
    u = (np.sqrt(rhoR) * uR + np.sqrt(rhoL) * uL) *1/(np.sqrt(rhoR) + np.sqrt(rhoL))
    v = (np.sqrt(rhoR) * vR + np.sqrt(rhoL) * vL) *1/(np.sqrt(rhoR) + np.sqrt(rhoL))
    h = (np.sqrt(rhoR) * hR + np.sqrt(rhoL) * hL) * 1/(np.sqrt(rhoR) + np.sqrt(rhoL))
    c = np.sqrt((gamma - 1)*(h - 0.5*(u**2 + v**2)))
    
    vn = u*etax + v*etay
    lambda1 = vn
    lambda2 = vn - c
    lambda3 = vn + c
    lambda4 = lambda1
    
    Rv = np.array([1, 1, 1, 0;
                   u - c*etax, u, u + a*etax, etay; 
                   v - c*etay, v, v + c*etay, -1*etax; 
                   h - c*vn, 0.5*(u**2 + v**2), h + c*vn, u*etay - v*etax])
    
def state_interpolation():
    pass


