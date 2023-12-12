import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii

direc = '/Users/jonty/mydata/stirred_discs/hd16743/scattered_light/vlt_sphere_julien/'

dstar = 57.932 #pc - convert arcsec to au for mass sensitivity

data = ascii.read(direc+'detection_limits.csv',guess=False,delimiter=',',data_start=2)

radius     = dstar*data["separation (arcsec)"].data
mass_limit = data["mass detection limit (Mjup)"].data
ml_low     = data["mass lower limit (Mjup)"].data
ml_hi      = data["mass upper limit (Mjup)"].data

rlimit = 300.

mask = np.where(radius <= rlimit)

radius = radius[mask]
mass_limit = mass_limit[mask]
ml_low = ml_low[mask]
ml_hi  = ml_hi[mask]

data2 = ascii.read(direc+'detection_limits_R1.csv',guess=False,delimiter=',',data_start=2)

radius2     = dstar*data2["separation (arcsec)"].data

mass_limit2 = data2["mass detection limit for 57 Myr (Mjup)"].data
ml_low2     = data2["mass lower limit for 38 Myr (Mjup)"].data
ml_hi2      = data2["mass upper limit for 76 Myr (Mjup)"].data

mask = np.where(radius2 <= rlimit)

radius2 = radius2[mask]
mass_limit2 = mass_limit2[mask]
ml_low2 = ml_low2[mask]
ml_hi2  = ml_hi2[mask]

radius3     = dstar*data2["separation (arcsec)"].data

mass_limit3 = data2["mass detection limit for 100 Myr (Mjup)"].data
ml_low3     = data2["mass lower limit for 40 Myr (Mjup)"].data
ml_hi3      = data2["mass upper limit for 600 Myr (Mjup)"].data

mask = np.where(radius3 <= rlimit)

radius3 = radius3[mask]
mass_limit3 = mass_limit3[mask]
ml_low3 = ml_low3[mask]
ml_hi3  = ml_hi3[mask]

fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(1, 1, 1)

disc_x = [118.0,118.0,198.0,198.0,118.0]
disc_y = [0.0,20.0,20.0,0.0,0.0]
ax.fill(disc_x,disc_y, fill=True, color="grey",alpha=0.1)

error_x = np.append(radius,radius[::-1])
error_y = np.append(ml_low,ml_hi[::-1])

ax.fill(error_x,error_y,linestyle='-',color='dodgerblue',alpha=0.1)
ax.plot(radius,mass_limit,linestyle='-',color='blue',label=r'5-$\sigma$ mass limit ($M_{\rm Jup}$), $t_{age} = 18^{+6}_{-6}~$Myr')

error_x2 = np.append(radius2,radius2[::-1])
error_y2 = np.append(ml_low2,ml_hi2[::-1])

ax.fill(error_x2,error_y2,linestyle='-',color='orange',alpha=0.1)
ax.plot(radius2,mass_limit2,linestyle='--',color='darkorange',label=r'5-$\sigma$ mass limit ($M_{\rm Jup}$), $t_{age} = 57^{+19}_{-19}~$Myr')

error_x3 = np.append(radius2,radius2[::-1])
error_y3 = np.append(ml_low3,ml_hi3[::-1])

ax.fill(error_x3,error_y3,linestyle='-',color='lightgreen',alpha=0.1)
ax.plot(radius3,mass_limit3,linestyle=':',color='green',label=r'5-$\sigma$ mass limit ($M_{\rm Jup}$), $t_{age} = 100^{+400}_{-60}~$Myr')

#single stirring planet
ax.errorbar([40.0],[4.0],xerr=[36.0],yerr=[10.0],linestyle='',marker='o',color='darkgray',label='Single Stirring Planet')
#single sculpting planet
ax.errorbar([40.0],[1.0],xerr=[40.0],yerr=[0.7],linestyle='',marker='o',color='black',label='Single Sculpting Planet')

ax.legend()
ax.set_xlabel(r'Radial separation (au)',fontsize=16)
ax.set_ylabel(r'Mass limit ($M_{\rm Jup}$)',fontsize=16)
ax.tick_params(axis="x", labelsize=12)
ax.tick_params(axis="y", labelsize=12)
ax.set_yscale('linear')
ax.set_xscale('log')
ax.set_ylim(0.0,20.0)
ax.set_xlim(0.0,300.0)
plt.tight_layout()
fig.savefig(direc + "HD16743_mass_limit.pdf",dpi=200)
plt.close()

