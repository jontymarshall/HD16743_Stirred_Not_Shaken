import numpy as np
import matplotlib.pyplot as plt
import time
from astropy.io import ascii
from scipy import interpolate
from RT_Code import RTModel

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

#plot the sed
def make_sed(m): 

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    
    try:
        ax.loglog(m.sed_wave, m.sed_emit, color='red',linestyle='-')
        ax.loglog(m.sed_wave, m.sed_warm, color='orange',linestyle='-')
        #for ij in range(0,int(m.parameters['nring'])):
        #    ax.loglog(m.sed_wave,m.sed_ringe[ij,:],linestyle='-',color='orange',alpha=0.1)
    except:
        print("No continuum emission model.")
    
    try:
        ax.loglog(m.sed_wave, m.sed_scat, color='blue',linestyle='-')
        #for ij in range(0,int(m.parameters['nring'])):  
        #    ax.loglog(m.sed_wave,m.sed_rings[ij,:],linestyle='-',color='dodgerblue',alpha=0.1)   
    except:
        print("No scattered light model.")

    ax.loglog(m.sed_wave, m.sed_emit+m.sed_warm, color='gray',linestyle='--')
    ax.loglog(m.sed_wave, m.sed_star, color='gray',linestyle=':')
    
    ax.loglog(m.sed_wave, m.sed_star + m.sed_emit + m.sed_scat + m.sed_warm, color='black',linestyle='-')

    try:
        ax.errorbar(m.obs_wave,m.obs_flux+yhot,yerr=m.obs_uncs,marker='o',linestyle='',mec='black',mfc='white',color='black')
    except:
        print("RTModel object has no observations.")

    ax.errorbar(m.obs_wave[-1],m.obs_flux[-1],yerr=m.obs_uncs[-1],marker='o',linestyle='',mec='darkred',mfc='red',color='red')
    ax.errorbar(irs_wav,irs_flx*1e3,xerr=None,yerr=irs_unc*1e3,marker=None,color='dodgerblue',alpha=0.1)
    
    ax.set_xlabel(r'Wavelength ($\mu$m)')
    ax.set_ylabel(r'Flux density (mJy)')
    ax.set_xlim(m.parameters["lmin"],m.parameters["lmax"])
    if np.max(m.sed_star) > np.max((m.sed_emit + m.sed_scat)):
        ax.set_ylim(10**(int(np.log10(np.max((m.sed_emit + m.sed_scat))) - 4)),10**(int(np.log10(np.max(m.sed_star)) + 1)))
    else:
        ax.set_ylim(10**(int(np.log10(np.max(m.sed_star)) - 4)),10**(int(np.log10(np.max((m.sed_emit + m.sed_scat))) + 1))) 
        
    fig.savefig(m.parameters['directory']+m.parameters['prefix']+'_sed.pdf',dpi=200)
    plt.close(fig)

    m.figure = fig

#Read in observations
import json
direc = '/Users/jonty/mydata/stirred_discs/hd16743/sed_fitting/'
data = json.load(open(direc+'hd16743.json'))

dstar = 1./data["main_results"][0]["plx_arcsec"]
lstar = data["main_results"][0]["lstar"] 

obs_photosphere = np.array(data['star_spec']['fnujy'])
obs_lambda = np.array(data['star_spec']['wavelength'])

data_lam  = np.array(data["phot_wavelength"][0])
data_flx  = np.array(data["phot_fnujy"][0])
data_unc  = np.array(data["phot_e_fnujy"][0])

data_ign = np.array(data["phot_ignore"][0])

for i in range(0,len(data_ign)):
    if data_lam[i] < 1.0 and data_flx[i] < 1.0:
        data_ign[i] = True
    print(data_ign[i],data_lam[i],data_flx[i])

data_lam = data_lam[np.where(data_ign == False)]
data_flx = data_flx[np.where(data_ign == False)]
data_unc = data_unc[np.where(data_ign == False)]



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

#add IRS spectrum
data_lam = np.append(data_lam,irs_l_interp)#,[1270.])
data_flx = np.append(data_flx,irs_f_interp)#,[1.2369e-3])
data_unc = np.append(data_unc,irs_u_interp)#,[0.1295e-3])
#add SPIRE 350um
data_lam = np.append(data_lam,[350.])
data_flx = np.append(data_flx,[0.0380])
data_unc = np.append(data_unc,[0.0060])
#add ALMA Band 6
data_lam = np.append(data_lam,[1270.])
data_flx = np.append(data_flx,[1.1790e-3])
data_unc = np.append(data_unc,[0.1180e-3])

import copy

xobs = copy.copy(data_lam[np.where((data_lam >= 0.4) & ((data_flx/data_unc) >= 3.))])
yobs = copy.copy(data_flx[np.where((data_lam >= 0.4) & ((data_flx/data_unc) >= 3.))])*1e3
uobs = copy.copy(data_unc[np.where((data_lam >= 0.4) & ((data_flx/data_unc) >= 3.))])*1e3
#print(xobs,yobs,uobs)

#Remove inner warm component from disc - replace with proper filters!
xwarm = np.array(data["model_spectra"][2]["wavelength"])
ywarm = np.array(data["model_spectra"][2]["fnujy"])*1e3

#print(xwarm,ywarm)

f = interpolate.interp1d(xwarm,ywarm)

yhot = f(xobs)

#print(yhot)

yobs = yobs - yhot

#benchmarking with time
start = time.time()

model = RTModel()
RTModel.get_parameters(model,'HD16743_RTModel_Input_File.txt')
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
model.sed_warm = f(model.sed_wave)
make_sed(model)

print("Fractional luminosity warm:",np.trapz(model.sed_warm,x=c/(model.sed_wave*1e-6))/np.trapz(model.sed_star,x=c/(1e-6*model.sed_wave)))
print("Fractional luminosity cold:",np.trapz(model.sed_emit-model.sed_warm,x=c/(model.sed_wave*1e-6))/np.trapz(model.sed_star,x=c/(1e-6*model.sed_wave)))

end = time.time()
multi_time = end - start
print("SED calculations took: ",multi_time," seconds.")