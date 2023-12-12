#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 18:52:48 2021

@author: jonty
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy import interpolate

c = 299792458.0 # m/s

direc = '/Users/jonty/mydata/stirred_discs/hd16743/co_spectrum/'

#Read in continuum image
hdul = fits.open(direc+'member.uid___A001_X1468_X1db.HD16743_sci.spw17_21_23_25.cont.I.pbcor.fits')
_image   = np.nan_to_num(hdul[0].data[0,0,:,:])
_iheadr  = hdul[0].header

pixl = _iheadr['CDELT2']*3600 #arcsec

x = pixl*(np.arange(0,_image.shape[0]) - int(_image.shape[0]/2))
y = pixl*(np.arange(0,_image.shape[1]) - int(_image.shape[1]/2))

xx,yy = np.meshgrid(x,y)

#Offsets of ring in x,y,z axes
deltax = 0.0
deltay = 0.0

#Angles for transformation
posang = -168.6*(np.pi/180) # degrees

c0 = np.cos(posang)
c1 = np.cos(0.0)
c2 = np.cos(0.5*np.pi)
s0 = np.sin(posang)
s1 = np.sin(0.0)
s2 = np.sin(0.5*np.pi)

#Transformations
trans1 = -(s0*c1*c2 + c0*s2)*yy + c2*deltax -s2*deltay
trans2 = -(s0*c1*s2 - c0*c2)*yy + s2*deltax +c2*deltay

#New x,y,z axes
x3= (c0*c1*c2 - s0*s2)*xx + trans1
y3= (c0*c1*s2 + s0*c2)*xx + trans2

#radius values
rr = np.sqrt(x3**2 + y3**2)

twosigma = 45e-6 #Jy/beam sigma noise level in continuum image

mask = np.zeros(_image.shape)
mask[np.where((_image >= twosigma)&(rr < 6.0))] = 1.

backmask = np.zeros(_image.shape)
backmask[np.where((rr >= 6.0)&(rr < 8.0))] = 1.


plt.imshow(mask,origin='lower') 

#Read in spw covering CO (2-1) line at 230.54 GHz
hdul = fits.open(direc+'member.uid___A001_X1468_X1db.HD16743_sci.spw25.cube.I.pbcor.fits')
_spect  = np.nan_to_num(hdul[0].data[0,:,:,:])
_sheadr = hdul[0].header

_spec = np.zeros(_spect.shape)
_specnoise = np.zeros(_spect.shape)

nu0 = _sheadr['CRVAL3'] #lowest frequency in spectral window 
dnu = _sheadr['CDELT3'] #frequency resolution of spectral window
nch = _sheadr['NAXIS3'] #number of channels

bmaj = _sheadr['BMAJ']*3600.0
bmin = _sheadr['BMIN']*3600.0

beam = np.pi *  (bmaj*bmin) / (4.0*np.log(2))

pixl = _sheadr['CDELT2']*3600 #arcsec

#calculate vorb
G = 6.67e-11 # kg m-2 s-1
mstar = 1.58 # M_sol
msol  = 1.99e30 # kg
cosi  = np.cos(87*(np.pi/180.0)) #disc inclination factor
au = 1.496e11 # m 

voffset = 14550.0 - 13737.669950575 #m/s - Gaia DR3 radial velocity + conversion from LSR to Helio, 21.9 km/s - Moor 2006

vkepl = cosi * np.sqrt( (G*mstar*msol) / (rr*au) )

vkepl[np.where(x3 < 0)] = -1.*vkepl[np.where(x3 < 0)]
vkepl[112,112] = 0.0

freqshift = 1 + ((vkepl-voffset)/c) 

plt.imshow(freqshift,origin='lower')

#extract spectra
frequency = nu0 + np.arange(0,nch)*dnu

channelwidth = dnu*c/2.30109e+11/1000

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.set_xlabel("Velocity shift (km/s)")
ax.set_ylabel("Flux density (mJy)")

velocities = (1+ (np.arange(-100e3,101e3,1e3)/c))*230.54e9

spec= np.zeros(len(velocities))
specnoise = np.zeros(len(velocities))

for i in range(0,_image.shape[0]):
    for j in range(0,_image.shape[1]):
        
        if mask[i,j] > 0.:
            
            freq = frequency*freqshift[i,j]
            
            f = interpolate.interp1d(freq,_spect[:,i,j])
            
            specshift = f(velocities)
            ax.plot(velocities,1000.*specshift,marker='',linestyle='-',color='dodgerblue',alpha=0.1)
            #print(np.min(1000*specshift),np.max(1000*specshift),np.std(1000*_spect[:,i,j]))
            spec += specshift*1000*pixl**2/beam/channelwidth
            

noise = np.zeros(nch)
for k in range(nch):
    noise[k] = np.std(_spect[k,np.where(backmask == 1)])*1000

f = interpolate.interp1d(frequency,noise)
specnoise = f(velocities)

meannoise = np.mean(specnoise)

print(meannoise)

velocities = (velocities - 230.54e9)*c/230.54e9/1e3
vmin = -50 # km/s
vmax = 50 #km/s

xf = np.asarray([-5,-5,+5,+5,-5]) + voffset*1e-3
yf = [-10,10,10,-10,-10]

ax.set_xlim(vmin,vmax)
ax.set_ylim(-5,10)
ax.errorbar(velocities,spec,yerr=specnoise,marker='',linestyle='-',color='blue',ecolor='dodgerblue')
ax.plot([vmin,vmax],[-2*meannoise,-2*meannoise],marker='',linestyle='--',color='darkgrey')
ax.plot([vmin,vmax],[2*meannoise,2*meannoise],marker='',linestyle='--',color='darkgrey')
ax.plot([vmin,vmax],[3*meannoise,3*meannoise],marker='',linestyle='--',color='darkgrey')
ax.plot([vmin,vmax],[4*meannoise,4*meannoise],marker='',linestyle='--',color='darkgrey')
ax.plot([vmin,vmax],[5*meannoise,5*meannoise],marker='',linestyle='--',color='darkgrey')
ax.fill(xf,yf,facecolor='lightblue',alpha=0.5)
ax.plot([vmin,vmax],[0,0],marker='',linestyle='-',color='black')
fig.savefig(direc+'HD16743_co_spectrum_R2.pdf',dpi=200)
plt.show()

#velocities = (np.arange(0,nch) - 0.5*nch)*channelwidth

#spec = np.zeros(nch)
#specnoise = np.zeros(nch)

#for each spectral channel
#for k in range(0,nch):
#    spec[k] = np.sum(_spect[k,np.where(mask == 1)])*1e3
#    specnoise[k] = np.std(_spect[k,np.where(backmask == 1)])*1e3
#    specbkgnd = np.mean(_spect[k,np.where(backmask == 1)])*1e3
#    spec[k] -= specbkgnd*(pixl**2*np.sum(mask)/beam)

#velocities = (velocities - 230.54e9)*c/230.54e9/1e3
#vmin = -50 # km/s
#vmax = 50 #km/s

#fig = plt.figure()
#ax = fig.add_subplot(1, 1, 1)
#ax.set_xlabel("Velocity shift (km/s)")
#ax.set_ylabel("Flux density (mJy/beam/km/s)")

#ax.set_xlim(vmin,vmax)
#ax.set_ylim(-5*np.median(specnoise),10*np.median(specnoise))
#ax.set_ylim(-5000,10000)
#ax.errorbar(velocities,spec,xerr=None,yerr=specnoise,marker='',linestyle='-',color='blue',ecolor='dodgerblue')
#ax.plot([vmin,vmax],[-1*specnoise,-1*specnoise],marker='',linestyle=':',color='lightgrey')
#ax.plot([vmin,vmax],[-2*specnoise,-2*specnoise],marker='',linestyle=':',color='lightgrey')
#ax.plot([vmin,vmax],[specnoise,specnoise],marker='',linestyle=':',color='darkgrey')
#ax.plot([vmin,vmax],[2*specnoise,2*specnoise],marker='',linestyle=':',color='darkgrey')
#ax.plot([vmin,vmax],[3*specnoise,3*specnoise],marker='',linestyle=':',color='darkgrey')
#ax.plot([vmin,vmax],[0,0],marker='',linestyle='-',color='black')
#fig.savefig(direc+'hd16743_co_spectrum.png',dpi=200)
#plt.show()
