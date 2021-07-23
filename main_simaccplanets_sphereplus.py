    from matplotlib import pyplot as plt
import numpy as np
import noise_model_hlines as nmh

#simulating a 1h observation of HIP65426
appmagJstar=6.826 #J mag of the star
JH_col=0.0 #J-H color of the star
distance=1e3/9.1566 #distance, parsec
DIT=60. #data integration time, second
NDIT=60. #number of integration
EWHa=10. #equivalent width of the emission line in Angstroms


appmagHstar=appmagJstar-JH_col

absmags=np.arange(12.,30.,0.1)

appmags=absmags+nmh.dist_modulus(distance)
contrastsJ=appmags-appmagJstar
contrastsH=appmags-appmagHstar

### Figure exploring the impact of the spectrograph type (lenslet vs fiber)
snrs_J_fb=np.zeros(len(absmags))
snrs_H_fb=np.zeros(len(absmags))
snrs_J_ls=np.zeros(len(absmags))
snrs_H_ls=np.zeros(len(absmags))

for k in range(0,len(appmags)):
    snrs_J_fb[k],_, _, _, _, _, _=nmh.noise_budget('J', appmags[k], contrastsJ[k], EWHa, 1.25, 50., DIT, NDIT, 283., 7., mode="fiber", nchannels=35*1000/35)
    snrs_H_fb[k],_, _, _, _, _, _=nmh.noise_budget('H', appmags[k], contrastsH[k], EWHa, 1.25, 50., DIT, NDIT, 283., 7., mode="fiber", nchannels=35*1000/35)
    snrs_J_ls[k],_, _, _, _, _, _=nmh.noise_budget('J', appmags[k], contrastsJ[k], EWHa, 1.25, 50., DIT, NDIT, 283., 7., mode="lenslet", nchannels=35*1000/35)
    snrs_H_ls[k],_, _, _, _, _, _=nmh.noise_budget('H', appmags[k], contrastsH[k], EWHa, 1.25, 50., DIT, NDIT, 283., 7., mode="lenslet", nchannels=35*1000/35)


fig, ax= plt.subplots()
ax.plot(absmags,snrs_J_fb,label='J, fiber')
ax.plot(absmags,snrs_H_fb,label='H, fiber')
ax.plot(absmags,snrs_J_ls,label='J, lenslets',linestyle='-')
ax.plot(absmags,snrs_H_ls,label='H, lenslets',linestyle='-')
ax.set_xlabel('Absolute magnitude')
ax.set_ylabel('S/N')
ax.set_ylim(0,50)
# ax.set_xlim=6
ax.axhline(y=5.,color='k',linestyle='--',label='5$\sigma$ threshold')
ax.set_title("Accreting planets around HIP65426, EW$_{Line}$=10$\AA$")
ax.legend()
fig.savefig('Fig_absmag_vs_SNR.pdf')


### Figure exploring the impact of the spectrograph spectral resolution
snrs_J_ls=np.zeros(len(absmags))
snrs_J_ls_1000=np.zeros(len(absmags))
snrs_J_ls_5000=np.zeros(len(absmags))
snrs_H_ls=np.zeros(len(absmags))
snrs_H_ls_1000=np.zeros(len(absmags))
snrs_H_ls_5000=np.zeros(len(absmags))

for k in range(0,len(appmags)):
    #J band
    snrs_J_ls[k],_, _, _, _, _, _=nmh.noise_budget('J', appmags[k], contrastsJ[k], EWHa, 1.25, 50., DIT, NDIT, 283., 7., mode="lenslet", nchannels=35.)
    snrs_J_ls_1000[k],_, _, _, _, _, _=nmh.noise_budget('J', appmags[k], contrastsJ[k], EWHa, 1.25, 1000., DIT, NDIT, 283., 7., mode="lenslet", nchannels=35.*1000/35)
    snrs_J_ls_5000[k],_, _, _, _, _, _=nmh.noise_budget('J', appmags[k], contrastsJ[k], EWHa, 1.25, 5000., DIT, NDIT, 283., 7., mode="lenslet", nchannels=35.*5000/35)

    #H band
    snrs_H_ls[k],_, _, _, _, _, _=nmh.noise_budget('H', appmags[k], contrastsH[k], EWHa, 1.25, 50., DIT, NDIT, 283., 7., mode="lenslet", nchannels=35.)
    snrs_H_ls_1000[k],_, _, _, _, _, _=nmh.noise_budget('H', appmags[k], contrastsH[k], EWHa, 1.25, 1000., DIT, NDIT, 283., 7., mode="lenslet", nchannels=35.*1000/35)
    snrs_H_ls_5000[k],_, _, _, _, _, _=nmh.noise_budget('H', appmags[k], contrastsH[k], EWHa, 1.25, 5000., DIT, NDIT, 283., 7., mode="lenslet", nchannels=35.*5000/35)


fig, ax= plt.subplots()
ax.plot(absmags,snrs_J_ls,label='J, R=50',color='b')
ax.plot(absmags,snrs_J_ls_1000,label='J, R=1000',linestyle='dashdot',color='b')
ax.plot(absmags,snrs_J_ls_5000,label='J, R=5000',linestyle=':',color='b')

ax.plot(absmags,snrs_H_ls,label='H, R=50',color='r')
ax.plot(absmags,snrs_H_ls_1000,label='H, R=1000',linestyle='dashdot',color='r')
ax.plot(absmags,snrs_H_ls_5000,label='H, R=5000',linestyle=':',color='r')

ax.set_xlabel('Absolute magnitude')
ax.set_ylabel('S/N')
ax.set_ylim(0,50)
ax.axhline(y=5.,color='k',linestyle='--',label='5$\sigma$ threshold')
ax.legend()
ax.set_title("Accreting planets around HIP65426, EW$_{Line}$=10$\AA$")
fig.savefig('Fig_absmag_vs_SNR_vs_R.pdf')
