# -*- coding: utf-8 -*-
"""
Created on Tue May  3 20:28:07 2022

@author: Dhruva
"""
import numpy as np

def Roe_Jacobian(curv, stateL, stateR, gamma = 1.4):
    #Want to operate using a Roe Averages Scheme of the Values
    #hat{A}
    
    etax = curv[0]
    etay = curv[1]
    
    rhoL = stateL[0]
    uL = stateL[1]/rhoL
    vL = stateL[2]/rhoL
    hL = (gamma-1)*(stateL[3] - 0.5*(uL**2 + vL**2))
    
    rhoR = stateR[0]
    uR = stateR[1]/rhoR
    vR = stateR[2]/rhoR
    hR = (gamma-1)*(stateR[3] - 0.5*(uR**2 + vR**2))
    
    rho = np.sqrt(rhoL * rhoR)
    u = (np.sqrt(rhoR) * uR + np.sqrt(rhoL) * uL) /(np.sqrt(rhoR) + np.sqrt(rhoL))
    v = (np.sqrt(rhoR) * vR + np.sqrt(rhoL) * vL) /(np.sqrt(rhoR) + np.sqrt(rhoL))
    h = (np.sqrt(rhoR) * hR + np.sqrt(rhoL) * hL) /(np.sqrt(rhoR) + np.sqrt(rhoL))
    c = np.sqrt(-(gamma - 1.)*(h - 0.5*(u**2 + v**2)))
    
    vn = u*etax + v*etay
    
    Rv = np.array([[1, 1, 1, 0],
                   [u - c*etax, u, u + c*etax, etay],
                   [v - c*etay, v, v + c*etay, -etax],
                   [h - c*vn, 0.5*(u**2 + v**2), h + c*vn, u*etay - v*etax]])
    print(Rv)
    eigval = [abs(vn - c), abs(vn), abs(vn + c), abs(vn)]
    Diag = np.zeros([len(eigval), len(eigval)])
    for i in range(len(eigval)):
        Diag[i, i] = eigval[i]
        
    g = (gamma-1)
    ss = c**2
    ek = 0.5*(u**2 + v**2)
    Lv = np.array([[(g*ek + c*vn)/(2.*ss), (-g*u - c*etax)*0.5/(ss), (-g*v - c*etay)*0.5/(ss), 0.5*g/(ss)],
                   [(ss - g*ek)/ss, g*u/ss, g*v/ss, -g/ss],
                   [(g*ek - c*vn)/(2.*ss), (-g*u + c*etax)/(2.*ss), (-g*v + c*etay)/(2.*ss), (-g/(2.*ss))],
                   [(v - vn*etay)/etax, etay, -etax, 0.]])
    R = np.matmul(Diag, Rv)
    out = np.matmul(Lv, R)
    print(out)
    return out
    
def Convective_Operator(curv, state, gamma=1.4):
    
    etax = curv[0]
    etay = curv[1]
    
    q1 = state[0]
    q2 = state[1]
    q3 = state[2]
    q4 = state[3]
    
    e1 = q2*etax + q3*etay
    e2 = (q2**2)/q1 *etax + q2*q3*etay/q1 + (gamma-1)*(q4 - 0.5*((q2**2)/q1 + (q3**2)/q1))*etax
    e3 = q2*q3*etax/q1 + (q3**2)*etay/q1 + (gamma-1)*(q4 - 0.5*((q2**2)/q1 + (q3**2)/q1))*etay
    e4 = gamma*q4 - 0.5*(gamma - 1)*((q2**2)/q1 + (q3**2)/q1)*(q2/q1*etax + q3/q1*etay)
    
    E = np.array([[e1], [e2], [e3], [e4]])
    
    return E

def phi(mode):
    if mode == "const":
        r = 1.
        phi = 1.
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
        Qnegone = state[i - 2, j]
        Qnegtwo = state[i - 4, j]
        Qzero = state[i + 2, j]
        Qposone = state[i + 4, j]
        
    if orient[direction] == 1:
        Qnegone = state[i, j - 2]
        Qnegtwo = state[i, j - 4]
        Qzero = state[i, j + 2]
        Qposone = state[i, j + 4]
        
    limiterL, rL = phi("const")
    QL = Qnegone + 0.25*epsilon*((1 - kappa)*(Qnegone - Qnegtwo) * limiterL + (1 + kappa)*(Q - Qnegone)/rL * limiterL)
    
    limiterR, rR = phi("const")
    QR = Qzero - 0.25*epsilon*((1 + kappa)*(Qzero - Qnegone)*limiterR*1/rR + (1 - kappa)*(Qposone - Qzero)*limiterR)
    
    return QL, QR
    

