#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 14:26:42 2021

@author: jonty
"""

from astropy.io import ascii
from scipy import interpolate
import numpy as np
import json

c = 2.99792458e8 # m/s

direc = '/Users/jonty/mydata/stirred_discs/hd16743/sed_fitting/'

#read in json file with photometry, stellar model, etc.
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
converters = {'wavelength': [ascii.convert_numpy(np.float32)],
              'flux': [ascii.convert_numpy(np.float32)],
              'error': [ascii.convert_numpy(np.float32)]}
irs_data = ascii.read(direc+'hd16743_spitzer_irs_spec_rescaled.tbl',converters=converters,format='csv',fast_reader=False)
irs_wav = irs_data["wavelength"].data
irs_flx = irs_data["flux"].data
irs_unc = irs_data["uncertainty"].data
irs_valid = np.where(irs_flx/irs_unc > 3.)

#append IRS spectro-photometry and ALMA Band 6 measurements to photometry
irs_l_interp = np.asarray([8.,12.,16.,20.,22.,24.,28.,32.,36.])
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
data_flx = np.append(data_flx,[1.2369e-3])
data_unc = np.append(data_unc,[0.1295e-3])

for i in range(0,len(data_lam)):
    print(data_lam[i],data_flx[i],data_unc[i])