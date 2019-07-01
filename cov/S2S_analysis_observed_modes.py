# coding: utf-8
"""This code takes summertime weekly prec data and analyzes the eof"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import mpl_toolkits.basemap as bm
from scipy import stats
import eofdata

# #### Observational data:####
# 
# GPCP (daily precipitation)
# 
# Extended austral summer: Nov-Mar
# Hindcast period: 1999/2000 to 2009/2010 (11 complete summers)
#Plot observed modes explained variance and lagged corrlation between modes
covariance = True
NWEEKS = 242
NX = 41
NY = 51
NMODES = 10
# Open observed data
PATH='/home/osman/datos/s2s_predictable_modes_paper/datos_S2S_Felipe/observation.gpcp/out/'
RUTA_OUT = '/home/osman/datos/s2s_predictable_modes_figures/cov/'
lon = np.arange(270, 270 + NX * 1.5, 1.5)
lat = np.arange(-60, -60 + NY * 1.5, 1.5)
[dx, dy] = np.meshgrid(lon, lat)
clevs = np.linspace(-40, 40, 11)
barra = plt.cm.bwr
for i in np.arange(1,5):
    WEEK = 'w' + str(i)
    file='weekly.anoconv.out.NOVMAR.' + WEEK +'.SA.bin'
    obs_weekly_precip = np.fromfile(PATH + file, dtype='f4', count=-1)
    obs_weekly_precip = np.reshape(obs_weekly_precip,[NWEEKS, NY * NX])
    #compute observed PC
    [total_var, obs_eval, obs_evec, obs_pc] = eofdata.eofdata(np.transpose(obs_weekly_precip),
                                                              NMODES, covariance)
    #plot first 4 modes
    fig = plt.figure(figsize=(22, 6), dpi=500)
    for j in np.arange(4):
        plt.subplot(1, 4, j + 1)
        mapproj = bm.Basemap(projection='cyl', llcrnrlon=270, llcrnrlat=-60,
                             urcrnrlon=270 + (NX - 1) *1.5,
                             urcrnrlat=-60 + (NY - 1) * 1.5, resolution='i')
        mapproj.drawcoastlines()
        lonproj, latproj = mapproj(dx, dy)
        CS1 = mapproj.contourf(lonproj, latproj, np.reshape(obs_evec[:, j], [NY, NX]), 
                               clevs, cmap=barra, vmin=-40,
                               vmax=40)
        plt.title("Mode " + str(j + 1), fontsize=12)
    #titulo general
    plt.suptitle('Observed Empirical Orthogonal Functions Week '+ str(i), fontsize=14, x=0.52,
                 y=0.95)
    cbar_ax = fig.add_axes([0.29, 0.05, 0.45, 0.05])
    fig.colorbar(CS1, cax=cbar_ax, orientation='horizontal')
    #cbar_ax.set_xticklabels([-0.9, -0.6, -0.3, 0, 0.3, 0.6, 0.9])#,size = 9)
    plt.savefig(RUTA_OUT + 'obs_eofs_week' + str(i) + '.png', dpi=500,
                bbox_inches='tight', orientation='landscape',
                papertype='A4') 

    #compute lagged correlation between PCs and plot
    lags = np.arange(0, 11, 1)
    signif = stats.norm.ppf(.975) / np.sqrt(242 - lags)
    fig = plt.figure(figsize = (10,17),dpi = 500)  #fig size in inches
    ax = plt.subplot(3,1,1)
    #lagged correlation btw PC1 and PC2
    ax.xcorr(obs_pc[0,:], obs_pc[1, :] , normed=True, usevlines=True, maxlags=10, linewidth=3)
    ax.axhline(y=0,xmin=0,linestyle='--',linewidth=1.2, color='b', alpha=0.6)
    ax.plot(lags, signif, linestyle='--', linewidth=1.2, color='b', alpha=0.6)
    ax.plot(lags, -signif, linestyle='--', linewidth=1.2, color='b', alpha=0.6)
    # tidy up the figure
    ax.set_xlim((0,11))
    ax.set_ylim((-0.5,0.5))
    plt.title('Lagged correlation btw PC1 and PC2',fontsize=10)
    #lagged correlation btw PC1 and PC3
    ax = plt.subplot(3,1,2)
    ax.xcorr(obs_pc[0,:], obs_pc[2, :] , normed=True, usevlines=True, maxlags=10, linewidth=3)
    ax.axhline(y=0,xmin=0,linestyle='--',linewidth=1.2, color='b', alpha=0.6)
    ax.plot(lags, signif, linestyle='--', linewidth=1.2, color='b', alpha=0.6)
    ax.plot(lags, -signif, linestyle='--', linewidth=1.2, color='b', alpha=0.6)
    # tidy up the figure
    ax.set_xlim((0,11))
    ax.set_ylim((-0.5,0.5))
    plt.title('Lagged correlation btw PC1 and PC3',fontsize=10)
    #lagged correlation btw PC1 and PC4
    ax = plt.subplot(3,1,3)
    ax.xcorr(obs_pc[0,:], obs_pc[3, :] , normed=True, usevlines=True, maxlags=10, linewidth=3)
    ax.axhline(y=0,xmin=0, linestyle='--', linewidth=1.2, color='b', alpha=0.6)
    ax.plot(lags, signif, linestyle='--', linewidth=1.2, color='b', alpha=0.6)
    ax.plot(lags, -signif, linestyle='--', linewidth=1.2, color='b', alpha=0.6)
    # tidy up the figure
    ax.set_xlim((0,11))
    ax.set_ylim((-0.5,0.5))
    plt.title('Lagged correlation btw PC1 and PC4',fontsize=10)
    fig.savefig(RUTA_OUT + 'PC1_laggedcorrelation_week' + str(i) + '.png',
                dpi=300, bbox_inches='tight', papertype='A4', orientation='landscape')

    #calculo la correlacion laggeada entre las Pcs 2-4
    fig = plt.figure(figsize = (10,17),dpi = 300)  #fig size in inches
    ax = plt.subplot(3,1,1)
    #lagged correlation btw PC2 and PC3
    ax.xcorr(obs_pc[1,:], obs_pc[2, :] , normed=True, usevlines=True, maxlags=10, linewidth=3)
    ax.axhline(y=0,xmin=0,linestyle='--',linewidth=1.2, color='b', alpha=0.6)
    ax.plot(lags, signif, linestyle='--', linewidth=1.2, color='b', alpha=0.6)
    ax.plot(lags, -signif, linestyle='--', linewidth=1.2, color='b', alpha=0.6)
    # tidy up the figure
    ax.set_xlim((0, 11))
    ax.set_ylim((-0.5, 0.5))
    plt.title('Lagged correlation btw PC2 and PC3', fontsize=10)
    ax = plt.subplot(3, 1, 2)
    #lagged correlation btw PC2 and PC4
    ax.xcorr(obs_pc[1, :], obs_pc[3, :] , normed=True, usevlines=True, maxlags=10, linewidth=3)
    ax.axhline(y=0, xmin=0, linestyle='--', linewidth=1.2, color='b', alpha=0.6)
    ax.plot(lags, signif, linestyle='--', linewidth=1.2, color='b', alpha=0.6)
    ax.plot(lags, -signif, linestyle='--', linewidth=1.2, color='b', alpha=0.6)
    # tidy up the figure
    ax.set_xlim((0, 11))
    ax.set_ylim((-0.5, 0.5))
    plt.title('Lagged correlation btw PC2 and PC4', fontsize=10)
    #lagged correlation btw PC3 and PC4
    ax = plt.subplot(3, 1, 3)
    ax.xcorr(obs_pc[2, :], obs_pc[3, :] , normed=True, usevlines=True, maxlags=10, linewidth=3)
    ax.axhline(y=0, xmin=0, linestyle='--', linewidth=1.2, color='b', alpha=0.6)
    ax.plot(lags, signif, linestyle='--', linewidth=1.2, color='b', alpha=0.6)
    ax.plot(lags, -signif, linestyle='--', linewidth=1.2, color='b', alpha=0.6)
    # tidy up the figure
    ax.set_xlim((0, 11))
    ax.set_ylim((-0.5, 0.5))
    plt.title('Lagged correlation btw PC3 and PC4', fontsize=10)
    fig.savefig(RUTA_OUT + 'PC2_laggedcorrelation_week' + str (i) +
                '.png', dpi=300, bbox_inches='tight', papertype='A4', orientation='landscape')


