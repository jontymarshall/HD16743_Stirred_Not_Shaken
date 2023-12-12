#!/usr/bin/env python3
#!/usr/bin/env bash
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 21 10:14:30 2018

@author: jonty
"""
import emcee
import corner
import dill
import json
import numpy as np
import matplotlib.pyplot as plt
import miepython.miepython as mpy
from pathos.multiprocessing import Pool
from numba import jit
import time
from astropy.io import ascii
from scipy import interpolate
from scipy import integrate

from RT_Code import RTModel

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 16}

plt.rc('font', **font)

direc = '/home/jmarshall/hd16743/'


#constants
h = 6.626e-34
c = 299792458.0 # m/s
k = 1.38e-23
sb = 5.67e-8 # 
au     = 1.495978707e11 # m 
pc     = 3.0857e16 # m
lsol   = 3.828e26 # W
rsol   = 6.96342e8 # m
MEarth = 5.97237e24 # kg

um = 1e-6 #for wavelengths in microns

#Plotting function for SED
def make_sed(m): 

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    
    try:
        ax.loglog(m.sed_wave, m.sed_emit, color='red',linestyle=':')
        #for ij in range(0,int(m.parameters['nring'])):
        #    ax.loglog(m.sed_wave,m.sed_ringe[ij,:],linestyle='-',color='orange',alpha=0.1)
    except:
        print("No continuum emission model.")
    
    try:
        ax.loglog(m.sed_wave, m.sed_scat, color='blue',linestyle=':')
        #for ij in range(0,int(m.parameters['nring'])):  
        #    ax.loglog(m.sed_wave,m.sed_rings[ij,:],linestyle='-',color='dodgerblue',alpha=0.1)   
    except:
        print("No scattered light model.")

    ax.loglog(m.sed_wave, (m.sed_emit + m.sed_scat), color='black',linestyle='--')
    ax.loglog(m.sed_wave, m.sed_star, color='black',linestyle='-.')
    
    ax.loglog(m.sed_wave, m.sed_star + m.sed_emit + m.sed_scat, color='black',linestyle='-')

    try:
        ax.errorbar(m.obs_wave,m.obs_flux,yerr=m.obs_uncs,marker='o',linestyle='',mec='black',mfc='white',color='black')
    except:
        print("RTModel object has no observations.")
    
    ax.set_xlabel(r'$\lambda$ ($\mu$m)')
    ax.set_ylabel(r'Flux density (mJy)')
    ax.set_xlim(m.parameters["lmin"],m.parameters["lmax"])
    if np.max(m.sed_star) > np.max((m.sed_emit + m.sed_scat)):
        ax.set_ylim(10**(np.log10(np.max((m.sed_emit + m.sed_scat))) - 6),10**(np.log10(np.max(m.sed_star)) + 1))
    else:
        ax.set_ylim(10**(np.log10(np.max(m.sed_star)) - 6),10**(np.log10(np.max((m.sed_emit + m.sed_scat))) + 1))  
    
    fig.savefig(m.parameters['directory']+m.parameters['prefix']+'_sed.png',dpi=200)
    plt.close(fig)

    m.figure = fig

#Read in data from json file
data = json.load(open(direc+'hd16743.json'))

dstar = 1./data["main_results"][0]["plx_arcsec"]
lstar = data["main_results"][0]["lstar"] 

obs_photosphere = np.array(data['star_spec']['fnujy'])
obs_lambda = np.array(data['star_spec']['wavelength'])


data_lam  = np.array(data["phot_wavelength"][0])
data_flx  = np.array(data["phot_fnujy"][0])
data_unc  = np.array(data["phot_e_fnujy"][0])

#convert photosphere for use in Hyperion
nu  = np.array(c/(obs_lambda*1e-4))
fnu = np.array(obs_photosphere)

#read in IRS spectrum
converters = {'wavelength': [ascii.convert_numpy(np.float64)],\
              'flux'      : [ascii.convert_numpy(np.float64)],\
              'error'     : [ascii.convert_numpy(np.float64)]}
irs_data = ascii.read(direc+'hd16743_spitzer_irs_spec_rescaled.tbl',format='csv',fast_reader=False)
irs_wav = irs_data["wavelength"].data
irs_flx = irs_data["flux"].data
irs_unc = irs_data["uncertainty"].data
irs_valid = np.where(irs_flx/irs_unc > 3.)

#append IRS spectro-photometry and ALMA Band 6 measurements to photometry
irs_l_interp = np.asarray([8.,12.,16.,22.,24.,28.,32.,36.])
f = interpolate.interp1d(irs_wav,irs_flx)
irs_f_interp = f(irs_l_interp)
f = interpolate.interp1d(irs_wav,irs_unc)
irs_u_interp = f(irs_l_interp)

#add IRS
data_lam = np.append(data_lam,irs_l_interp)#,[1270.])
data_flx = np.append(data_flx,irs_f_interp)#,[1.2369e-3])
data_unc = np.append(data_unc,irs_u_interp)#,[0.1295e-3])
#add SPIRE
data_lam = np.append(data_lam,[250.,350.])
data_flx = np.append(data_flx,[0.0820,0.0380])
data_unc = np.append(data_unc,[0.0060,0.0060])
#add ALMA
data_lam = np.append(data_lam,[1270.])
data_flx = np.append(data_flx,[1.1790e-3])*1e3
data_unc = np.append(data_unc,[0.1180e-3])*1e3

import copy

xobs = copy.copy(data_lam[np.where((data_lam >= 0.4) & ((data_flx/data_unc) >= 3.))])
yobs = copy.copy(data_flx[np.where((data_lam >= 0.4) & ((data_flx/data_unc) >= 3.))])
uobs = copy.copy(data_unc[np.where((data_lam >= 0.4) & ((data_flx/data_unc) >= 3.))])
#print(xobs,yobs,uobs)

#Remove inner warm component from disc - replace with proper filters!
xwarm = np.array(data["model_spectra"][2]["wavelength"])
ywarm = np.array(data["model_spectra"][2]["fnujy"])

#print(xwarm,ywarm)

f = interpolate.interp1d(xwarm,ywarm*1e3)
yhot = f(xobs)

#print(yhot)

yobs = yobs - yhot

#print("Photometry:")
#print("-----------")
#for i in range(0,len(xobs)):
#    print(xobs[i],yobs[i],uobs[i])

#emcee depends on you having defined a reasonable likelihood and prior, 
#everything else is just brute force. We will assume flat priors today, but I 
#suggest you read the relevant sections of Hogg+2010 to understand why this 
#assumption should be discarded whenever possible.
def lnprior(theta):
    psdex,grains,mdust=theta
    
    if -4.0 < psdex < -3.0 and -0.5 < grains < 1.5 and -6 < mdust < 3:
        return 0.0
    return -np.inf

def lnlike(theta,xobs,yobs,uobs):
    psdex,grains,mdust=theta

    model = RTModel()
    RTModel.get_parameters(model,'HD16743_RTModel_Input_File.txt')
    
    #Apply dust values to the model
    model.parameters['mdust'] = 10**mdust
    model.parameters['q'] = psdex
    model.parameters['amin'] = 10**grains
    model.obs_wave = xobs
    model.obs_flux = yobs
    model.obs_uncs = uobs
    RTModel.make_star(model)
    RTModel.scale_star(model,lrange=[0.6,9.0])
    RTModel.make_dust(model)
    RTModel.make_disc(model)
    RTModel.calculate_surface_density(model)
    RTModel.read_optical_constants(model)
    RTModel.calculate_qabs(model)
    RTModel.calculate_dust_emission(model,mode='full',tolerance=0.02)
    RTModel.calculate_dust_scatter(model)
    RTModel.flam_to_fnu(model)

    #make_sed(model,(-1,-1,-1))
    
    f = interpolate.interp1d(model.sed_wave,(model.sed_star+model.sed_emit))
    sed_model = f(model.obs_wave) 
    

    el = np.where(model.obs_wave >= 10.)
    #print(psdex,10**grains,10**mdust,np.sum(((model.obs_flux[el] - sed_model[el])**2/model.obs_uncs[el]**2)))

    return -0.5 * np.sum(((model.obs_flux[el] - sed_model[el])**2/model.obs_uncs[el]**2))

    
def lnprob(theta,xobs,yobs,uobs):
    lp=lnprior(theta)
    
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, xobs,yobs,uobs)

def run_emcee(sampler,ndim,labels,steps=100,prefix="",state=None):
    print("Running MCMC...")
    state = sampler.run_mcmc(state,steps,rstate0=np.random.get_state(),progress=True)
    print("Done.")
    
    plt.clf()
    fig, axes = plt.subplots(ndim, 1, sharex=True, figsize=(8, 9))
    
    
    for i in range(ndim):
        axes[i].plot(sampler.chain[:, :, i].T, color="k", alpha=0.4)
        axes[i].set_ylabel(labels[i])
    
    fig.tight_layout(h_pad=0.0)
    fig.savefig(prefix+"line-time.png")
    return sampler,state

def mcmc_results(sampler,ndim,percentiles=[16, 50, 84],burnin=200,labels="",prefix=""):
    
    samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))
    print(samples.shape)
    
    fig = corner.corner(samples, labels=labels[0:ndim])
    fig.savefig(prefix+"line-triangle.png")
    credible_interval=[]
    for i in range(ndim):
        credible_interval.append(np.percentile(samples[:,i], percentiles))
        credible_interval[i][2] -= credible_interval[i][1]
        credible_interval[i][0] = credible_interval[i][1] - credible_interval[i][0]
    
    print("MCMC results:")
    for i in range(ndim):
        print("{0}  = {1[1]} + {1[2]} - {1[0]}".format(labels[i],credible_interval[i]))
    
    print("Run RT model for maximum probability SED:")    
    #small dust grains

    dustmass = credible_interval[2][1] #M_Earth
    dustsize = credible_interval[1][1] #um
    dustqval = credible_interval[0][1] #
    
    #Set up SED modelling
    model = RTModel()

    RTModel.get_parameters(model,'HD16743_RTModel_Input_File.txt')

    model.parameters['mdust'] = 10**dustmass
    model.parameters['q']     = dustqval
    model.parameters['amin']  = 10**dustsize
    
    model.obs_wave = xobs
    model.obs_flux = yobs
    model.obs_uncs = uobs

    RTModel.make_star(model)
    RTModel.scale_star(model,lrange=[0.6,9.0])
    RTModel.make_dust(model)
    RTModel.make_disc(model)
    RTModel.calculate_surface_density(model)
    RTModel.read_optical_constants(model)
    RTModel.calculate_qabs(model)
    RTModel.calculate_dust_emission(model,mode='full',tolerance=0.02)
    RTModel.calculate_dust_scatter(model)
    RTModel.flam_to_fnu(model)

    make_sed(model)
    
    print("Calculated model chi-squared...")
    f = interpolate.interp1d(model.sed_wave,(model.sed_star+model.sed_emit))
    sed_model = f(model.obs_wave) 

    el = np.where(model.obs_wave >= 10.)
    #print(psdex,10**grains,10**mdust,np.sum(((model.obs_flux[el] - sed_model[el])**2/model.obs_uncs[el]**2)))

    csq = np.sum(((model.obs_flux[el] - sed_model[el])**2/model.obs_uncs[el]**2))
     
    print("Model chi-squared: ",csq)
    
    print("Written model SED output...")
    ascii.write([model.sed_wave,model.sed_star,model.sed_disc,(model.sed_star+model.sed_disc)], prefix+'_astrosil_sed_component.tbl',\
                 names=['wavelength','fstar','fdust','ftotal'], overwrite=True)
    
    print("Plot maximum probability model SED...")
    make_sed(model)
    print("Done.")

#SET UP AND RUN EMCEE
nwalkers = 30
nsteps = 300
nburn = int(0.8*nsteps)
ndim = 3
#radius,rwidth,psdex,grains,rho0=theta
pos=[[np.random.uniform(-3.8,-3.2),np.random.uniform(-0.5,0.5),np.random.uniform(-2,0)] for w in range(nwalkers)]
labels=[r"$\gamma$",r"$a_{\rm min}$",r"$M_{\rm dust}$"]


#Instantiate RTModel object
starttime = time.time()

sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, pool=Pool(80), args=(xobs,yobs,uobs))
sampler,state = run_emcee(sampler,ndim,labels,nsteps,prefix=direc+"HD16743_MCMC_SED_Fitting_",state=pos)
dill.dump(sampler,open(direc+"HD16743_MCMC_SED_Fitting_Sampler_Output_nw"+str(int(nwalkers))+"_ns"+str(int(nsteps))+".p","wb"))

mcmc_results(sampler,ndim,labels=labels,prefix=direc+"HD16743_MCMC_SED_Fitting_nw"+str(int(nwalkers))+"_ns"+str(int(nsteps))+"_",burnin=nburn)

endtime = time.time()
print("Processing took ",endtime-starttime," seconds.")