# -*- coding: utf-8 -*-
"""
Created on Mon May  2 21:44:40 2022

@author: Dhruva
"""

from GridReader import reader, metrics, halo_augmenter
import matplotlib.pyplot as plt

file = "Fine.dat"
data, grid, haloXL, haloXR, haloYT, haloYB = reader(file)
gridh = halo_augmenter(grid, haloXL, haloXR, haloYB, haloYT)
total = metrics(gridh)


fig = plt.figure(dpi=500)
ax1 = fig.add_subplot(111)
ax1.scatter(data[:, 0], data[:, 1], s = 0.1)
ax1.scatter(haloXL[0], haloXL[1],  s=0.3, c="g")
ax1.scatter(haloXR[0], haloXR[1], s=0.3, c="g")
ax1.scatter(haloYT[0], haloYT[1], s=0.3, c="r")
ax1.scatter(haloYB[0], haloYB[1], s=0.3, c="r")
ax1.set_xlabel(r"$\xi$")
ax1.set_ylabel(r"$\eta$")

fig2, ax2 = plt.subplots(dpi=500)
ax2.set_aspect('equal')
tcf = ax2.tricontourf(total[1::2, 1::2, 0].flatten(), total[1::2, 1::2, 1].flatten(), total[1::2, 1::2, 2].flatten(), 25)
fig2.colorbar(tcf)
# ax2.scatter(center[:, 0], center[:, 1], volume)
ax2.set_xlabel(r"$\xi$")
ax2.set_ylabel(r"$\eta$")

fig3, ax3 = plt.subplots(dpi=500)
ax3.set_aspect('equal')
#Eta fluxes
tcf = ax3.tricontourf(total[0::2, 1::2, 0].flatten(), total[0::2, 1::2, 1].flatten(), total[0::2, 1::2, 2].flatten(), 500)
fig3.colorbar(tcf)
# ax3.title("$\eta$ Fluxes")
ax3.set_xlabel(r"$\xi$")
ax3.set_ylabel(r"$\eta$")

fig4, ax4 = plt.subplots(dpi=500)
# ax3.set_aspect('equal')
#Xi fluxes
tcf = ax4.tricontourf(total[1::2, 0::2, 0].flatten(), total[1::2, 0::2, 1].flatten(), total[1::2, 0::2, 2].flatten(), 500)
fig4.colorbar(tcf)
# ax4.title("$\eta$ Fluxes")
ax4.set_xlabel(r"$\xi$")
ax4.set_ylabel(r"$\eta$")