#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  26 16:25:48 2022

@author: jonty
"""

import numpy as np
import matplotlib.pyplot as plt


mstar = 1.537 # M_sol
rdisc = 157.7   # au
wdisc = 79.7    # au
inc   = 0.13  # radians
tsys  = 18    # Myrs
tsys_s = tsys*1e6*365.25*24*3600 #s
mearth = 6e24 # kg
msun   = 2e30 # kg
au     = 1.496e11 # m
G      = 6.67e-11
rho    = 1000 #kg/m3

sigma_ = np.logspace(-3,4,num=50,endpoint=True,base=10.0) # kg / m2
mass_  = np.logspace(-6,3,num=50,endpoint=True,base=10.0) # M_Earth

vrel = 103.9*mstar**0.5*inc*rdisc**(-0.5)*1e3 # m/s
vesc = (G * (mstar*msun) / (rdisc*au) )**0.5  # m/s

Adisc = np.pi*au**2*((rdisc+0.5*wdisc)**2 - (rdisc-0.5*wdisc)**2) #au2

#surface density limits taken from Matra et al. (2019)
sigma_eqn12 = 4.43e3 * mstar**1.5 * rdisc**(-0.5) * inc**4 / tsys / mass_

rbody = ((3*mass_*mearth)/(4*np.pi*rho))**(1/3)
sig = 4.*np.pi*rbody**2
nbody = (1./tsys_s) / sig / vrel / (1. + (vesc**2/vrel**2))

sigma_eqn13 = (mass_*mearth) * nbody

print(nbody,sigma_eqn13)

drdisc = 8 * 3**0.5 * (rdisc*au) * ((mass_*mearth) / (mstar*msun))**(1/3)
sigma_eqn14 = (mass_*mearth) / (2 * np.pi * (rdisc*au) * drdisc)

#mass limits
#vesc = vrel
mbody = mass_*mearth
rbody = (0.24 * mbody / rho)**(1/3)
vbody = (G * mbody / rbody)**(1/2)
delta = (vbody - vrel)**2

mass_vels = mbody[np.argmin(delta)]/mearth

#clearing



#plot
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(1, 1, 1)
ax.set_xlabel(r'Mass of largest body ($M_{\oplus}$)',fontsize=16)
ax.set_ylabel(r'Belt surface density (kg m$^{-2}$)',fontsize=16)
ax.set_ylim(1e-3,1e4)
ax.set_xlim(1e-6,1e3)
ax.loglog(mass_,sigma_eqn12,marker='',linestyle='-',color='purple')
#ax.loglog(mass_,sigma_eqn13,marker='v',linestyle='',color='orange')
ax.loglog(mass_,sigma_eqn14,marker='^',linestyle='',color='green')
ax.loglog([mass_vels,mass_vels],[sigma_[0],sigma_[-1]],marker='',linestyle='-.',color='blue')
ax.loglog([mass_[0],mass_[-1]],[0.3,0.3],marker='',linestyle=':',color='red')
plt.show()
