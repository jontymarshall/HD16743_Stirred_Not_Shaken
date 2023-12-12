#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 19:13:38 2022

@author: jonty
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits


direc = '/Users/jonty/mydata/stirred_discs/hd16743/scattered_light/hst_nicmos_elodie/'

files = ['f110w/HIP-12361_NICMOS_F110W_MRDI_planet-mode_Lib-281_KL-93_Combined_smoothed.fits',
         'f110w_5p/Best_Forward_Model_F110W.fits',
         'f160w/HIP-12361_NICMOS_F160W_MRDI_planet-mode_Lib-397_KL-132_Combined_smoothed.fits',
         'f160w_5p/Best_Forward_Model_F160W.fits']

#F110W
hdul = fits.open(direc+files[0])
hdr = hdul[0].header
img = hdul[0].data
img = np.nan_to_num(img,nan=99,posinf=99,neginf=99)

nx = int(hdr['NAXIS2']) # Image size keyworks from FITS header
ny = int(hdr['NAXIS1']) # Image size keyworks from FITS header

img = img[6:nx-6,6:ny-6]

nx -= 12
ny -= 12

pixscale = hdr['PIXSIZE']
photfnu = hdr['PHOTFNU']

img *= 1e6 * photfnu / pixscale**2 #in Jy/arcsec**2

rpeak = 153.
rfwhm =  88.

rmin = rpeak - (rfwhm/2)
rmax = rpeak + (rfwhm/2)

pa  = 79.*(np.pi/180.)
inc = 82.*(np.pi/180.)
dstar = 57.932
#calculate model grid for dust density calculation 
xlabel = ['+5','+4','+3','+2','+1','0','-1','-2','-3','-4','-5']
ylabel = ['-5','-4','-3','-2','-1','0','+1','+2','+3','+4','+5']
nlabel = len(xlabel)
xindex = ((nx-1)/(nlabel-1))*np.arange(0,11)
yindex = ((ny-1)/(nlabel-1))*np.arange(0,11)

fig = plt.figure(figsize=(6,7))
ax = fig.add_subplot(1, 1, 1)
ax.set_ylabel(r'North Offset ($^{\prime\prime}$)',fontsize=16)
ax.set_xlabel(r'East Offset ($^{\prime\prime}$)',fontsize=16)
ax.tick_params(axis="x", labelsize=12)
ax.tick_params(axis="y", labelsize=12)
plt.xticks(xindex,xlabel,fontsize=12)
plt.yticks(yindex,ylabel,fontsize=12)
im = ax.imshow(img,vmin=-10,vmax=80,origin='lower',cmap='gnuplot2')

cbaxes = fig.add_axes([0.13, 0.85, 0.77, 0.03])
cb = plt.colorbar(mappable=im,cax=cbaxes,orientation="horizontal",ticklocation='top')
cb.set_label(r'Flux ($\mu$Jy/arcsec$^{2}$)', size=18)
cb.ax.tick_params(labelsize=12)

fig.savefig(direc+'HD16743_hst-nicmos_f110W_image.pdf',dpi=200)
plt.show()
plt.close()

hdul = fits.open(direc+files[1])
hdr = hdul[0].header
mod = hdul[0].data
mod = np.nan_to_num(mod,nan=99,posinf=99,neginf=99)

nx = int(hdr['NAXIS2']) # Image size keyworks from FITS header
ny = int(hdr['NAXIS1']) # Image size keyworks from FITS header
try:
    pixscale = hdr['PIXSIZE']
except: 
    pixscale = 0.07565150000000001
try:
    photfnu = hdr['PHOTFNU']
except:
    photfnu = 1.21121E-06

mod = mod[6:nx-6,6:ny-6]

nx -= 12
ny -= 12

mod *= 1e6 * photfnu / pixscale**2 #in Jy/arcsec**2

fig = plt.figure(figsize=(6,7))
ax = fig.add_subplot(1, 1, 1)
ax.set_ylabel(r'North Offset ($^{\prime\prime}$)',fontsize=16)
ax.set_xlabel(r'East Offset ($^{\prime\prime}$)',fontsize=16)
ax.tick_params(axis="x", labelsize=12)
ax.tick_params(axis="y", labelsize=12)
plt.xticks(xindex,xlabel,fontsize=12)
plt.yticks(yindex,ylabel,fontsize=12)
im = ax.imshow(mod,vmin=-10,vmax=80,origin='lower',cmap='gnuplot2')

cbaxes = fig.add_axes([0.13, 0.85, 0.77, 0.03])
cb = plt.colorbar(mappable=im,cax=cbaxes,orientation="horizontal",ticklocation='top')
cb.set_label(r'Flux ($\mu$Jy/arcsec$^{2}$)', size=18)
cb.ax.tick_params(labelsize=12)

fig.savefig(direc+'HD16743_hst-nicmos_f110W_model.pdf',dpi=200)
plt.show()
plt.close()

res = np.zeros(img.size)
res = img - mod
mask = np.where(img >= 99)
res[mask] = mod[mask]


fig = plt.figure(figsize=(6,7))
ax = fig.add_subplot(1, 1, 1)
ax.set_ylabel(r'North Offset ($^{\prime\prime}$)',fontsize=16)
ax.set_xlabel(r'East Offset ($^{\prime\prime}$)',fontsize=16)
ax.tick_params(axis="x", labelsize=12)
ax.tick_params(axis="y", labelsize=12)
plt.xticks(xindex,xlabel,fontsize=12)
plt.yticks(yindex,ylabel,fontsize=12)
im = ax.imshow(res,vmin=-10,vmax=80,origin='lower',cmap='gnuplot2')

cbaxes = fig.add_axes([0.13, 0.85, 0.77, 0.03])
cb = plt.colorbar(mappable=im,cax=cbaxes,orientation="horizontal",ticklocation='top')
cb.set_label(r'Flux ($\mu$Jy/arcsec$^{2}$)', size=18)
cb.ax.tick_params(labelsize=12)

fig.savefig(direc+'HD16743_hst-nicmos_f110W_residual.pdf',dpi=200)
plt.show()
plt.close()

#F160W
hdul = fits.open(direc+files[2])
hdr = hdul[0].header
img = hdul[0].data
img = np.nan_to_num(img,nan=99,posinf=99,neginf=99)

nx = int(hdr['NAXIS2']) # Image size keyworks from FITS header
ny = int(hdr['NAXIS1']) # Image size keyworks from FITS header

img = np.pad(img,2,'constant',constant_values=99)

nx += 4
ny += 4

pixscale = hdr['PIXSIZE']
photfnu = hdr['PHOTFNU']

img *= 1e6 * photfnu / pixscale**2 #in Jy/arcsec**2

rpeak = 153.
rfwhm =  88.

rmin = rpeak - (rfwhm/2)
rmax = rpeak + (rfwhm/2)

pa  = 79.*(np.pi/180.)
inc = 82.*(np.pi/180.)
dstar = 57.932
#calculate model grid for dust density calculation 
xlabel = ['+5','+4','+3','+2','+1','0','-1','-2','-3','-4','-5']
ylabel = ['-5','-4','-3','-2','-1','0','+1','+2','+3','+4','+5']
nlabel = len(xlabel)
xindex = ((nx-1)/(nlabel-1))*np.arange(0,11)
yindex = ((ny-1)/(nlabel-1))*np.arange(0,11)

fig = plt.figure(figsize=(6,7))
ax = fig.add_subplot(1, 1, 1)
ax.set_ylabel(r'North Offset ($^{\prime\prime}$)',fontsize=16)
ax.set_xlabel(r'East Offset ($^{\prime\prime}$)',fontsize=16)
ax.tick_params(axis="x", labelsize=12)
ax.tick_params(axis="y", labelsize=12)
plt.xticks(xindex,xlabel,fontsize=12)
plt.yticks(yindex,ylabel,fontsize=12)
im = ax.imshow(img,vmin=-10,vmax=80,origin='lower',cmap='gnuplot2')

cbaxes = fig.add_axes([0.13, 0.85, 0.77, 0.03])
cb = plt.colorbar(mappable=im,cax=cbaxes,orientation="horizontal",ticklocation='top')
cb.set_label(r'Flux ($\mu$Jy/arcsec$^{2}$)', size=18)
cb.ax.tick_params(labelsize=12)

fig.savefig(direc+'HD16743_hst-nicmos_f160W_image.pdf',dpi=200)
plt.show()
plt.close()

hdul = fits.open(direc+files[3])
hdr = hdul[0].header
mod = hdul[0].data
mod = np.nan_to_num(mod,nan=99,posinf=99,neginf=99)

nx = int(hdr['NAXIS2']) # Image size keyworks from FITS header
ny = int(hdr['NAXIS1']) # Image size keyworks from FITS header
try:
    pixscale = hdr['PIXSIZE']
except:
    pixscale = 0.07565150000000001

try:
    photfnu = hdr['PHOTFNU']
except:
    photfnu = 1.49585E-06

mod = np.pad(mod,2,mode='constant',constant_values=99)
nx += 4
ny += 4

mod *= 1e6 * photfnu / pixscale**2 #in Jy/arcsec**2

fig = plt.figure(figsize=(6,7))
ax = fig.add_subplot(1, 1, 1)
ax.set_ylabel(r'North Offset ($^{\prime\prime}$)',fontsize=16)
ax.set_xlabel(r'East Offset ($^{\prime\prime}$)',fontsize=16)
ax.tick_params(axis="x", labelsize=12)
ax.tick_params(axis="y", labelsize=12)
plt.xticks(xindex,xlabel,fontsize=12)
plt.yticks(yindex,ylabel,fontsize=12)
im = ax.imshow(mod,vmin=-10,vmax=80,origin='lower',cmap='gnuplot2')

cbaxes = fig.add_axes([0.13, 0.85, 0.77, 0.03])
cb = plt.colorbar(mappable=im,cax=cbaxes,orientation="horizontal",ticklocation='top')
cb.set_label(r'Flux ($\mu$Jy/arcsec$^{2}$)', size=18)
cb.ax.tick_params(labelsize=12)

fig.savefig(direc+'HD16743_hst-nicmos_f160W_model.pdf',dpi=200)
plt.show()
plt.close()

res = np.zeros(img.size)
res = img - mod
mask = np.where(img >= 99)
res[mask] = mod[mask]

fig = plt.figure(figsize=(6,7))
ax = fig.add_subplot(1, 1, 1)
ax.set_ylabel(r'North Offset ($^{\prime\prime}$)',fontsize=16)
ax.set_xlabel(r'East Offset ($^{\prime\prime}$)',fontsize=16)
ax.tick_params(axis="x", labelsize=12)
ax.tick_params(axis="y", labelsize=12)
plt.xticks(xindex,xlabel,fontsize=12)
plt.yticks(yindex,ylabel,fontsize=12)
im = ax.imshow(res,vmin=-10,vmax=80,origin='lower',cmap='gnuplot2')

cbaxes = fig.add_axes([0.13, 0.85, 0.77, 0.03])
cb = plt.colorbar(mappable=im,cax=cbaxes,orientation="horizontal",ticklocation='top')
cb.set_label(r'Flux ($\mu$Jy/arcsec$^{2}$)', size=18)
cb.ax.tick_params(labelsize=12)

fig.savefig(direc+'HD16743_hst-nicmos_f160W_residual.pdf',dpi=200)
plt.show()
plt.close()