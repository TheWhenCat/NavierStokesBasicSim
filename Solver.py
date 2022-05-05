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
    
    Rv = np.array([[1, 1, 1, 0],
                   [u - c*etax, u, u + c*etax, etay],
                   [v - c*etay, v, v + c*etay, -1*etax],
                   [h - c*vn, 0.5*(u**2 + v**2), h + c*vn, u*etay - v*etax]])
    
    eigval = [abs(vn - c), abs(vn), abs(vn + c), abs(vn)]
    Diag = np.zero([len(eigval), len(eigval)])
    for i in len(range(eigval)):
        Diag[i, i] = eigval[i]
        
    g = (gamma-1)
    ss = c**2
    ek = 0.5*(u**2 + v**2)
    Lv = np.array([[(g*ek + c*vn)/(2*ss), (-1*g*u - c*etax)*0.5*1/(ss), (-1*g*v - c*etay)*0.5*1/(ss), 0.5*g/(ss)],
                   [(ss - g*ek)/ss, g*u/ss, g*v/ss, -1*g/ss],
                   [(g*ek - c*vn)/(2*ss), (-1*g*u + c*etax)/(2*ss), (-1*g*v + c*etay)/(2*ss), (-1*g/(2*ss))],
                   [(v - vn*etay)/etax, etay, -etax, 0]])
    R = np.matmul(Diag, Lv)
    out = np.matmul(Lv, Diag)
    return out
    
def Convective_OperatorE():
    pass

def Convective_OperatorF():
    pass

def phi(mode):
    if mode == "const":
        r = 1
        phi = 1
        return phi,r
    if mode == "minmod":
        return AttributeError("Dumbass You Need to Fix This")
    
    
def state_interpolation(index, epsilon, kappa, direction, state):
    orient = {"xi":0, "eta":1}
    #i index is flipped in my grid
    i = index[0]
    j = index[1]
    Q = state[i, j]
    if orient[direction] == 0:
        Qnegone = state[i - 1, j]
        Qnegtwo = state[i - 3, j]
        Qposone = state[i + 1, j]
        
    limiterL, rL = phi("const")
    QL = Qnegone + 0.25*epsilon*((1 - kappa)*(Qnegone - Qnegtwo) * limiterL + (1 + kappa)*(Q - Qnegone)*1/rL * limiterL)
    
    limiterR, rR = phi("const")
    QR = Q - 0.25*epsilon*((1 + kappa)*(Q - Qnegone)*limiterR*1/rR + (1 - kappa)*(Qposone - Q)*limiterR)
    return QL, QR
    

