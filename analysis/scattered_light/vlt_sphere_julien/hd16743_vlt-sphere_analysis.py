#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 11:43:22 2020

@author: jonty
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.ndimage.interpolation import rotate 

direc = '/Users/jonty/mydata/stirred_discs/hd16743/scattered_light/'

#read in fits file
#H-band image
#orientation N up, E left
#flux units Jy/arcsec^2
#plate scale 12.5 mas per pixel
hdul = fits.open(direc+'vlt_sphere_julien/HD16743_799x799_sum_sadi_vip_filtered_80_spiders_smoothed_Jy_per_arcsec2.fits')
hdr = hdul[0].header
img = hdul[0].data
img = np.nan_to_num(img) 
nx = int(hdr['NAXIS2']) # Image size keyworks from FITS header
ny = int(hdr['NAXIS1']) # Image size keyworks from FITS header

rpeak = 153.
rfwhm =  88.

rmin = rpeak - (rfwhm/2)
rmax = rpeak + (rfwhm/2)

pa  = 79.*(np.pi/180.)
inc = 82.*(np.pi/180.)
dstar = 57.932

pixscale = 12.5e-3 # 12.5 mas per pixel

#calculate model grid for dust density calculation 
rscale = pixscale*dstar
nxc = int(nx/2)
nyc = int(ny/2)

imgmask = np.zeros(img.shape)

x =  rscale*(np.arange(0,nx)  - nxc)
y =  rscale*(np.arange(0,ny)  - nyc)

xx,yy = np.meshgrid(x,y)
rr = np.sqrt(xx**2 + yy**2)

xp = xx*np.cos(pa) - yy*np.sin(pa)
yp = yy*np.cos(pa) + xx*np.sin(pa)
rvmax = (xp**2/rmax**2) + (yp**2/(rmax*np.cos(inc))**2) 
imgmask[np.where(rvmax <= 1.)] = 1.
rvmin = (xp**2/rmin**2) + (yp**2/(rmin*np.cos(inc))**2 )
imgmask[np.where(rvmin < 1.)] = 0.
imgmask[np.where(rr <= 125)] = 0.

annulus = np.where(imgmask > 0)

print(np.sum(img[annulus]),np.median(img[annulus]),np.std(img[annulus]))

rotang = np.arange(0,19*np.pi/18.,step=np.pi/18)

theta = []
surfb = []

for ang in rotang: 
    imgmask = np.zeros(img.shape)

    x =  rscale*(np.arange(0,nx)  - nxc)
    y =  rscale*(np.arange(0,ny)  - nyc)
    
    xx,yy = np.meshgrid(x,y)
    rr = np.sqrt(xx**2 + yy**2)
    
    xp = xx*np.cos(ang) - yy*np.sin(ang)
    yp = yy*np.cos(ang) + xx*np.sin(ang)
    rvmax = (xp**2/rmax**2) + (yp**2/(rmax*np.cos(inc))**2) 
    imgmask[np.where(rvmax <= 1.)] = 1.
    rvmin = (xp**2/rmin**2) + (yp**2/(rmin*np.cos(inc))**2 )
    imgmask[np.where(rvmin < 1.)] = 0.
    imgmask[np.where(rr <= 125)] = 0.
    
    annulus = np.where(imgmask > 0)
    
    #plt.imshow(imgmask)
    
    theta.append(ang)
    surfb.append(np.median(img[annulus]))
    
    #print(ang,np.sum(img[annulus]),np.median(img[annulus]),np.std(img[annulus]))

theta = np.asarray(theta)
surfb = np.asarray(surfb)*1e6

#plot angular surface brightness profile
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(1, 1, 1)
ax.set_ylabel(r'Median Flux ($\mu$Jy/arcsec$^{2}$)',fontsize=16)
ax.set_xlabel(r'Position Angle ($^{o}$)',fontsize=16)
ax.tick_params(axis="x", labelsize=12)
ax.tick_params(axis="y", labelsize=12)
ax.set_ylim(-5,10)
ax.set_yscale('linear')
ax.set_xlim(0,180)
ax.set_xscale('linear')
ax.plot(theta*(180./np.pi),surfb,linestyle='-',marker='o',mec='blue',mfc='dodgerblue',color='blue')
ax.plot([90.0,90.0],[-1e3,1e3],linestyle='--',color='gray')
#ax.plot([100.,100.],[-1e3,1e3],linestyle='--',color='red')
ax.plot([0,180],[0,0],linestyle='-',color='black')
xlabel = ['-90','-70','-50','-30','-10','+10','+30','+50','+70','+90']
nlabel = len(xlabel)
xindex = ((180)/(nlabel-1))*np.arange(0,nlabel)
plt.xticks(xindex,xlabel,fontsize=14)
plt.tight_layout()
fig.savefig(direc+'HD16743_SPHERE_angular_profile.pdf',dpi=200)
plt.show()
plt.close()

radial = np.arange(0.5*rfwhm,2.5*rpeak + 0.5*rfwhm,step=0.1*rfwhm)

radis = []
surfs = []
surfu = []

beamarea = np.pi * 52e-3**2 / 4. / np.log(2)

for rad in radial: 
    imgmask = np.zeros(img.shape)

    x =  rscale*(np.arange(0,nx)  - nxc)
    y =  rscale*(np.arange(0,ny)  - nyc)
    
    xx,yy = np.meshgrid(x,y)
    rr = np.sqrt(xx**2 + yy**2)
    
    rout = rad + (rfwhm/2)
    rin  = rad - (rfwhm/2)
    
    if rin < 0:
        rin = 0.
    
    xp = xx*np.cos(pa) - yy*np.sin(pa)
    yp = yy*np.cos(pa) + xx*np.sin(pa)
    rvmax = (xp**2/rout**2) + (yp**2/(rout*np.cos(inc))**2) 
    imgmask[np.where(rvmax <= 1.)] = 1.
    rvmin = (xp**2/rin**2) + (yp**2/(rin*np.cos(inc))**2 )
    imgmask[np.where(rvmin < 1.)] = 0.
    
    annulus = np.where(imgmask > 0)
    nonannulus = np.where((imgmask != 1)&(rr <= rout)&(rr > rin))
    #plt.imshow(imgmask)
    
    annarea = pixscale**2 * np.pi*(rout**2 - rin**2)
    #scalearea = np.sqrt(annarea/beamarea)
    radis.append(rad)
    surfs.append(np.sum(img[annulus])/annarea)
    surfu.append(5.0*np.std(img[nonannulus]))#/scalearea)
    print(rad*pixscale,np.sum(img[annulus])/annarea,np.std(img[annulus]))

radis = np.asarray(radis)
surfs = np.asarray(surfs)*1e6
surfu = np.asarray(surfu)
#plot angular surface brightness profile
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(1, 1, 1)
ax.set_ylabel(r'Median Flux ($\mu$Jy/arcsec$^{2}$)',fontsize=16)
ax.set_xlabel(r'Angular distance from star ($^{\prime\prime}$)',fontsize=16)
ax.tick_params(axis="x", labelsize=12)
ax.tick_params(axis="y", labelsize=12)
ax.set_ylim(0.0,100.0)
ax.set_yscale('linear')
ax.set_xlim(0.0,5.0)
ax.set_xscale('linear')
ax.plot(radis/dstar,surfs,linestyle='-',marker='o',mec='blue',mfc='dodgerblue',color='blue')
ax.plot(radis/dstar,surfu,linestyle='-',marker='',color='orange')
ax.plot([rpeak/dstar,rpeak/dstar],[-1000,1000],linestyle='-',color='darkgray')
ax.plot([(rpeak+0.5*rfwhm)/dstar,(rpeak+0.5*rfwhm)/dstar],[-1000,1000],linestyle=':',color='gray')
ax.plot([(rpeak-0.5*rfwhm)/dstar,(rpeak-0.5*rfwhm)/dstar],[-1000,1000],linestyle=':',color='gray')

ax.plot([-90,90],[0.0,0.0],linestyle='-',color='black')
plt.tight_layout()
fig.savefig(direc+'HD16743_SPHERE_radial_profile.pdf',dpi=200)
plt.show()
plt.close()


#plot image
hdul = fits.open(direc+'vlt_sphere_julien/HD16743_799x799_sum_sadi_vip_filtered_80_spiders_smoothed_Jy_per_arcsec2.fits')
hdr = hdul[0].header
img = hdul[0].data*1e6
img = np.nan_to_num(img)

fig, ax = plt.subplots(figsize=(8,10))
img_scaled = img
img_scaled[np.where(rr <100)] = 99
img_scaled[np.where(rr >250)] = 99
im = ax.imshow(img_scaled,vmin=-20,vmax=50,cmap='gnuplot2',origin='lower')
xsu = 798/4
ax.set_xticks([0,xsu,2*xsu,3*xsu,4*xsu])
ax.set_xticklabels(['+5.0','+2.5','0','-2.5','-5.0'],fontsize=18)
ax.set_yticks([0,xsu,2*xsu,3*xsu,4*xsu])
ax.set_yticklabels(['-5.0','-2.5','0','+2.5','+5.0'],fontsize=18)
ax.set_xlabel(r'East Offset (")',fontsize=24)
ax.set_ylabel(r'North Offset (")',fontsize=24)
#Add star
#ax.plot(45.,45.,"*",color="black",markersize=18)
#Add beam
beamx = 50. + 3.4*4.16*np.sin((180/np.pi)*np.arange(253))
beamy = 50. + 3.4*4.16*np.cos((180/np.pi)*np.arange(253))
ax.fill(beamx,beamy, fill=True, color="black",alpha=1.0)



#Add disc ellipse
from matplotlib.patches import Ellipse
router = rpeak+rfwhm/2
ellipse = Ellipse((400,400),2*router/(dstar*pixscale),2*router*np.cos(inc)/(dstar*pixscale),angle = 79.0,color='white',alpha=0.5)
ax.add_patch(ellipse)


#Add colour bar
cbaxes = fig.add_axes([0.15, 0.85, 0.81, 0.03])
cb = plt.colorbar(mappable=im,cax=cbaxes,cmap='gnuplot2',orientation="horizontal",ticklocation='top')
cb.set_label(r'Flux ($\mu$Jy/arcsec$^{2}$)', size=24)
cb.ax.tick_params(labelsize=18)
plt.tight_layout()
outputdirec = '/Users/jonty/mydata/stirred_discs/hd16743/scattered_light/'
fig.savefig(outputdirec + "HD16743_vlt-sphere_image_smoothed.pdf",dpi=200)



#plot artefacts
hdul = fits.open(direc+'vlt_sphere_julien/HD16743_799x799_synthetic_psf_residuals_derotated_stacked_filtered_80_spiders.fits')
hdr = hdul[0].header
art = hdul[0].data
art = np.nan_to_num(art) 
nx = int(hdr['NAXIS2']) # Image size keyworks from FITS header
ny = int(hdr['NAXIS1']) # Image size keyworks from FITS header

pixscale = 12.5e-3 # 12.5 mas per pixel

fig, ax = plt.subplots(figsize=(8,10))

art_scaled = art/pixscale**2
art_scaled[np.where(rr >250)] = 99
art_scaled[np.where(art_scaled < -5)] = -5
art_scaled[np.where(art_scaled > 5)] = 5
ar = ax.imshow(art_scaled,vmin=-2,vmax=5,cmap='gnuplot2',origin='lower')
xsu = 798/4
ax.set_xticks([0,xsu,2*xsu,3*xsu,4*xsu])
ax.set_xticklabels(['+5.0','+2.5','0','-2.5','-5.0'],fontsize=18)
ax.set_yticks([0,xsu,2*xsu,3*xsu,4*xsu])
ax.set_yticklabels(['-5.0','-2.5','0','+2.5','+5.0'],fontsize=18)
ax.set_xlabel(r'East Offset (")',fontsize=24)
ax.set_ylabel(r'North Offset (")',fontsize=24)
#Add star
#ax.plot(45.,45.,"*",color="black",markersize=18)
#Add beam
beamx = 50. + 4.16*np.sin((180/np.pi)*np.arange(253))
beamy = 50. + 4.16*np.cos((180/np.pi)*np.arange(253))
ax.fill(beamx,beamy, fill=True, color="white",alpha=0.5)


#Add colour bar
cbaxes = fig.add_axes([0.15, 0.85, 0.81, 0.03])
cb = plt.colorbar(mappable=ar,cax=cbaxes,cmap='gnuplot2',orientation="horizontal",ticklocation='top')
cb.set_label(r'Intensity (arb. units)', size=24)
cb.ax.tick_params(labelsize=18)
plt.tight_layout()
outputdirec = '/Users/jonty/mydata/stirred_discs/hd16743/scattered_light/'
fig.savefig(outputdirec + "HD16743_vlt-sphere_residuals.pdf",dpi=200)

#plot artefacts
hdul = fits.open(direc+'vlt_sphere_julien/HD16743_799x799_fake_disk_convolved.fits')
hdr = hdul[0].header
model = hdul[0].data
model = np.nan_to_num(model) 
nx = int(hdr['NAXIS2']) # Image size keyworks from FITS header
ny = int(hdr['NAXIS1']) # Image size keyworks from FITS header

pixscale = 12.5e-3 # 12.5 mas per pixel

fig, ax = plt.subplots(figsize=(8,10))
model_scaled = model/(0.25*np.max(model))
model_scaled[np.where(rr >250)] = 99

ar = ax.imshow(model_scaled,vmin=0.0,vmax=1.0,cmap='gnuplot2',origin='lower')
xsu = 798/4
ax.set_xticks([0,xsu,2*xsu,3*xsu,4*xsu])
ax.set_xticklabels(['+5.0','+2.5','0','-2.5','-5.0'],fontsize=18)
ax.set_yticks([0,xsu,2*xsu,3*xsu,4*xsu])
ax.set_yticklabels(['-5.0','-2.5','0','+2.5','+5.0'],fontsize=18)
ax.set_xlabel(r'East Offset (")',fontsize=24)
ax.set_ylabel(r'North Offset (")',fontsize=24)
#Add star
#ax.plot(45.,45.,"*",color="black",markersize=18)
#Add beam
beamx = 50. + 3.4*4.16*np.sin((180/np.pi)*np.arange(253))
beamy = 50. + 3.4*4.16*np.cos((180/np.pi)*np.arange(253))
ax.fill(beamx,beamy, fill=True, color="white",alpha=0.5)


#Add colour bar
cbaxes = fig.add_axes([0.15, 0.85, 0.81, 0.03])
cb = plt.colorbar(mappable=ar,cax=cbaxes,cmap='gnuplot2',orientation="horizontal",ticklocation='top')
cb.set_label(r'Intensity (arb. units)', size=24)
cb.ax.tick_params(labelsize=18)
plt.tight_layout()
outputdirec = '/Users/jonty/mydata/stirred_discs/hd16743/scattered_light/'
fig.savefig(outputdirec + "HD16743_vlt-sphere_model.pdf",dpi=200)