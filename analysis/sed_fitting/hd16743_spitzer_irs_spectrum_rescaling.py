#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 17:21:48 2020

@author: jonty
"""

import astropy
from astropy.io import ascii
from scipy import interpolate
import matplotlib.pyplot as plt
import numpy as np
import json

#read in json file with photometry, stellar model, etc.
direc = '/Users/jonty/mydata/stirred_discs/hd16743/sed_fitting/'
data = json.load(open(direc+'hd16743.json'))

dstar = 1./data["main_results"][0]["plx_arcsec"]
lstar = data["main_results"][0]["lstar"] 

obs_photosphere = np.array(data['star_spec']['fnujy'])
obs_lambda = np.array(data['star_spec']['wavelength'])

data_lam  = np.array(data["phot_wavelength"][0])
data_flx  = np.array(data["phot_fnujy"][0])
data_unc  = np.array(data["phot_e_fnujy"][0])

#read in IRS spectrum
converters = {'wavelength': [ascii.convert_numpy(np.float32)],
              'flux': [ascii.convert_numpy(np.float32)],
              'error': [ascii.convert_numpy(np.float32)]}
irs_data = ascii.read(direc+'cassis_tbl_spcf_15018240t.tbl',converters=converters,format='ipac',fast_reader=False)
irs_wav = irs_data["wavelength"].data
irs_flx = irs_data["flux"].data
irs_unc1 = irs_data["errrms"].data
irs_unc2 = irs_data["errsys"].data
irs_unc3 = irs_data["errcal"].data
irs_unc = (irs_unc1**2 + irs_unc2**2 + irs_unc3**2 )**0.5
irs_order = irs_data["strip"].data
irs_valid = np.where(irs_flx/irs_unc > 3.)

#plot up original IRS spectrum with orders in different colours
colours = ['red','orange','yellow','green','blue','purple']
irs_orders = np.unique(irs_order)

#scale IRS spectrum
f = interpolate.interp1d(obs_lambda,obs_photosphere)
star_flux = f(irs_wav)

csq_min = 1e30

for scalefactor in np.arange(1.0,1.2,0.0001):
    csq = np.sum((((scalefactor*irs_flx[np.where(irs_wav < 10.)]) - star_flux[np.where(irs_wav < 10.)])/irs_unc[np.where(irs_wav < 10.)])**2)
    
    if csq < csq_min:
        csq_min = csq
        scf_min = scalefactor

#plot up rescaled IRS spectrum vs. stellar photosphere and original
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(1, 1, 1)
ax.set_xlabel(r'Wavelength [$\mu$m]',fontsize=16)
ax.set_ylabel(r'Flux density [Jy]',fontsize=16)
ax.set_ylim(0.05,.30)
ax.set_xlim(5.,38.0)
ax.loglog(obs_lambda,obs_photosphere,color='black')
for order in irs_orders:
    value = np.where(irs_order == order)
    if order == 1:
        irs_flx[value] = scf_min*irs_flx[value]
    ax.errorbar(irs_wav[value],irs_flx[value],color=colours[order-1],xerr=None,yerr=irs_unc[value])

plt.show()

#irs_flx[np.where(irs_order == 1)] = scf_min*irs_flx[np.where(irs_order == 1)]


index = np.argsort(irs_wav)

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(1, 1, 1)
ax.set_xlabel(r'Wavelength [$\mu$m]',fontsize=16)
ax.set_ylabel(r'Flux density [Jy]',fontsize=16)
ax.set_ylim(0.05,.30)
ax.set_xlim(5.,38.0)
ax.errorbar(irs_wav[index],irs_flx[index],yerr=irs_unc[index],color='black',linestyle='-')
plt.show()


#save rescaled IRS spectrum 
from astropy.table import Table, Column

data = Table([irs_wav[index],irs_flx[index],irs_unc[index]], names=['wavelength', 'flux', 'uncertainty'])
ascii.write(data, 'hd16743_spitzer_irs_spec_rescaled.tbl', overwrite=True,delimiter=',')
