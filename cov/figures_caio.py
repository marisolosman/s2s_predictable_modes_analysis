# # Austral summer sub-seasonal precipitation predictable modes over South America #
# 
# ## Questions to be addressed: ##
# + How do the analyzed S2S models (and the multi-model ensemble of these models) represent the dominant modes of sub-seasonal precipitation variability during the extended austral summer at different lead times? 
# - What modes of sub-seasonal precipitation variability contribute to the potential predictability during the austral summer? 
# * What percentage of the actual S2S models prediction ability comes from the identified modes of sub-seasonal precipitation variability? 
# 

# In[2]:
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import mpl_toolkits.basemap as bm
from scipy import stats
covariance = True
WEEK = ['w1', 'w2', 'w3']
WEEK = ['w4']
RUTA = './outputs/'
RUTA_OUT = './figures_caio/'
NWEEKS = 242
NX = 41
NY = 51
NMODES = 3
NPOINTS = NX*NY
models = ['cma', 'eccc', 'ecmwf', 'ncep', 'mmem']
#np.savez('./outputs/EOF_analsys_w1.npz',obs_pc_w1=obs_pc,obs_evec_w1=obs_evec,
#         models_pc_w1=models_pc_sorted, models_evec_w1=models_evec_sorted, dx=dx, dy=dy)
#clevs = np.linspace(-40, 40, 21)
#minv = clevs[0]
#maxv = clevs[-1]
#barra = plt.cm.bwr #colorbar
##obs week 1
#FILE = RUTA + 'EOF_analsys_' + WEEK[0] + '.npz'
#f = np.load(FILE)
#var_obs = f['obs_evec_' + WEEK[0]]
#dx = f['dx']
#dy = f['dy']
#fig = plt.figure(figsize=(26, 10), dpi=300)
#for i in range(4):
#    plt.subplot(1, 4, 1 + i)
#    mapproj = bm.Basemap(projection='cyl', llcrnrlon=270, llcrnrlat=-60,
#                         urcrnrlon=270 + (NX - 1) *1.5, urcrnrlat=-60 + (NY - 1) * 1.5,
#                         resolution='i')
#    mapproj.drawcoastlines()
#    lonproj, latproj = mapproj(dx, dy)
#    var = np.reshape(var_obs[i, :], [NY, NX])
#    CS1 = mapproj.contourf(lonproj, latproj, var, clevs, cmap=barra, vmin=minv, vmax=maxv)
#    plt.title("Mode " + str(i + 1), fontsize=24)
#plt.suptitle('Empirical Orthogonal Functions', fontsize=30, x=0.51, y=0.92)
#plt.subplots_adjust(right=0.9, top=0.90)#, wspace=None, hspace=None)
#cbar_ax = fig.add_axes([0.92, 0.29, 0.03, 0.45])
#cax=fig.colorbar(CS1, cax=cbar_ax, orientation='vertical')
#ticklabs = cax.ax.get_yticklabels()
#cax.ax.set_yticklabels(ticklabs, fontsize=24)
#plt.savefig(RUTA_OUT + 'obs_eofs.png', dpi=300,
#            bbox_inches='tight', orientation='portrait', papertype='A4') 
#plt.close()
#f = np.load(RUTA + 'pss_exp_var_w1.npz')
#obs_var = f['arr_1']
#error_eval = np.sqrt(2 / NWEEKS) * obs_var * 100
#fig = plt.figure
#plt.errorbar(np.arange(1,11), obs_var * 100, error_eval,
#             color='b', linewidth=1.5)
#plt.xlabel('Modes', fontsize=14)
#plt.ylabel('percentage variance (%)', fontsize=14)
#plt.title('Scree plot - Observed precipitation', fontsize=20)
#plt.savefig(RUTA_OUT + 'variance_eofs.png', dpi=300, bbox_inches='tight',
#            orientation='landscape', papertype='A4')
#
#figure 3 to 5 EFO1 to 3 for Obs and w1-w4 forecast
#for i in range(3):
#    fig = plt.figure(figsize=(26, 42), dpi=300)
#    for j in range(4):
#        signo = -1 if np.logical_and(j == 2, i==0) else 1 
#        FILE = RUTA + 'EOF_analsys_' + WEEK[j] + '.npz'
#        f = np.load(FILE)
#        var_obs = f['obs_evec_' + WEEK[j]]
#        var_mod = f['models_evec_' + WEEK[j]]
#
#        dx = f['dx']
#        dy = f['dy']
#        plt.subplot(6, 4, 1 + j)
#        mapproj = bm.Basemap(projection='cyl', llcrnrlon=270, llcrnrlat=-60,
#                             urcrnrlon=270 + (NX - 1) *1.5, urcrnrlat=-60 + (NY - 1) * 1.5,
#                             resolution='i')
#        mapproj.drawcoastlines()
#        lonproj, latproj = mapproj(dx, dy)
#        var = signo * np.reshape(var_obs[i, :], [NY, NX])
#        CS1 = mapproj.contourf(lonproj, latproj, var, clevs, cmap=barra, vmin=minv, vmax=maxv)
#        plt.title("Obs", fontsize=24)
#        for k in range(5):
#            #print(models[k], np.corrcoef(var_mod[k, i, :], var_obs[i])[0,1])
#
#            plt.subplot(6, 4, 5 + j + k * 4)
#            mapproj = bm.Basemap(projection='cyl', llcrnrlon=270, llcrnrlat=-60,
#                                 urcrnrlon=270 + (NX - 1) *1.5, urcrnrlat=-60 + (NY - 1) * 1.5,
#                                 resolution='i')
#            mapproj.drawcoastlines()
#            lonproj, latproj = mapproj(dx, dy)
#            var =  signo * np.reshape(var_mod[k, i, :], [NY, NX])
#            CS1 = mapproj.contourf(lonproj, latproj, var, clevs, cmap=barra, vmin=minv, vmax=maxv)
#            plt.title(models[k].upper() + ' - Week ' + str(j + 1), fontsize=24)
##        titulo general
#    plt.suptitle('Empirical Orthogonal Function: '+ str(i + 1), fontsize=30, x=0.51, y=0.92)
#    plt.subplots_adjust(right=0.9, top=0.90)#, wspace=None, hspace=None)
#    cbar_ax = fig.add_axes([0.92, 0.33, 0.05, 0.35])
#    cax=fig.colorbar(CS1, cax=cbar_ax, orientation='vertical')
#    ticklabs = cax.ax.get_yticklabels()
#    cax.ax.set_yticklabels(ticklabs, fontsize=24)
#    plt.savefig(RUTA_OUT + 'sorted_eofs_' + str(i + 1) + '.png', dpi=300,
#                bbox_inches='tight', orientation='portrait', papertype='A4') 
#
#x = np.arange(NWEEKS)
#for i in range(3):
#    fig = plt.figure(figsize=(41, 27), dpi=300)  #fig size in inches
#    for j in np.arange(4):
#        signo = -1 if np.logical_and(j == 2, i==0) else 1 
#        FILE = RUTA + 'EOF_analsys_' + WEEK[j] + '.npz'
#        f = np.load(FILE)
#        color = iter(cm.rainbow(np.linspace(0, 1, 6)))
#        c = next(color)
#        ax = plt.subplot(4, 1, j+1)
#        y_obs = signo * f['obs_pc_' + WEEK[j]][i, :]
#        ax.plot(x, y_obs, color=c, linewidth=1.5, label='Obs')
#        for k in np.arange(5):
#            c = next(color)
#            y = signo * f['models_pc_' + WEEK[j]][k, i, :]
#            #print(models[k], str(j+1), np.corrcoef(y,y_obs)[0,1])
#            ax.plot(x, y, color=c, label=models[k].upper(), linewidth=0.8)
#        ax.axhline(y=0, xmin=0, linestyle='--', linewidth=0.8, color='b', alpha=0.5)
#        # tidy up the figure
#        ax.set_xlim((0, NWEEKS+1))
#        ax.set_ylim((-5, 5))
#        plt.tick_params(axis='both', which='major', labelsize=18)
#        plt.title('WEEK ' + str(j+1), fontsize=24)
#        #plt.legend(loc='lower right', fontsize=24)
#    plt.subplots_adjust(right=0.9, top=0.90)
#    plt.subplots_adjust(wspace=0.05) 
#    plt.legend(bbox_to_anchor=(1.02, 2), loc='lower left', borderaxespad=0., fontsize=24)
#    plt.suptitle('PC ' + str(i + 1), fontsize=28, x=0.52, y=0.93)
#    fig.savefig(RUTA_OUT + 'sorted_pcs_' + str(i + 1) +'.png', dpi=300, bbox_inches='tight',
#                papertype='A4', orientation='landscape')
#
##
#### Plotting results
### ### Potential predictability
##
### In[13]:
#clevs = np.linspace(-1, 1, 21)
#barra = plt.cm.bwr #colorbar
#for i in [0]:
#    f = np.load('./outputs/skill_predic_' + WEEK[i] + '_bis.npz')
#    R_forec_predic_modes = f['arr_0']
#    R_forec_residual_modes = f['arr_1']
#    R_forec_anom = f['arr_2']
#    R_obs_predic_modes = f['arr_3']
#    print(np.min(R_forec_predic_modes))
#    print(np.max(R_forec_predic_modes))
#    print(R_forec_predic_modes)
#
#    print('obs', np.mean(np.reshape(R_obs_predic_modes, NX*NY)))
#    #f['R_obs_predic_modes1']
#    dx = f['arr_5']
#    dy = f['arr_6']
#    fig = plt.figure(figsize=(30, 28), dpi=300)
#    plt.subplot(4, 5, 3)
#    mapproj = bm.Basemap(projection='cyl',llcrnrlon = 270, llcrnrlat = -60,
#                         urcrnrlon = 270+(NX-1)*1.5, urcrnrlat = -60+(NY-1)*1.5,
#                         resolution = 'i')
#    mapproj.drawcoastlines()
#    lonproj, latproj = mapproj(dx, dy)      #poject grid
#    CS1 = mapproj.contourf(lonproj, latproj, R_obs_predic_modes[:, :], clevs,
#                           cmap=barra,vmin=-1,vmax=1) #extended generate pretty colorbar
#    plt.title('Observed predict. modes', fontsize=22)
#    for k in range(5):
#        print(models[k], 'Anom', np.mean(np.reshape(R_forec_anom[k, :, :], NX*NY)))
#        print(models[k], ' PM ', np.mean(np.reshape(R_forec_predic_modes[k, :, :], NX*NY)))
#        print(models[k], ' Res ', np.mean(np.reshape(R_forec_residual_modes[k, :, :], NX*NY)))
#        plt.subplot(4, 5,k + 6)
#        mapproj = bm.Basemap(projection='cyl',llcrnrlon = 270, llcrnrlat = -60,
#                             urcrnrlon = 270+(NX-1)*1.5, urcrnrlat = -60+(NY-1)*1.5,
#                             resolution = 'i')
#        mapproj.drawcoastlines()
#        lonproj, latproj = mapproj(dx, dy)      #poject grid
#        CS1 = mapproj.contourf(lonproj, latproj, R_forec_anom[k, :, :], clevs,
#                               cmap=barra,vmin=-1,vmax=1) #extended generate pretty colorbar
#        plt.title('Forecasted anom. - ' + models[k].upper(), fontsize=22)
#        plt.subplot(4, 5, k + 11)
#        mapproj = bm.Basemap(projection='cyl', llcrnrlon=270, llcrnrlat=-60,
#                             urcrnrlon=270 + (NX - 1) * 1.5, urcrnrlat=-60 + (NY - 1) * 1.5,
#                             resolution='i')
#        mapproj.drawcoastlines()
#        lonproj, latproj = mapproj(dx, dy)      #poject grid
#        CS1 = mapproj.contourf(lonproj, latproj, R_forec_predic_modes[k, :, :], clevs,
#                               cmap=barra, vmin=-1, vmax=1) #extended generate pretty colorbar
#        plt.title('Forec. predict. modes - ' + models[k].upper(), fontsize=22)
#        plt.subplot(4, 5, k + 16)
#        mapproj = bm.Basemap(projection='cyl', llcrnrlon=270, llcrnrlat=-60,
#                             urcrnrlon=270 + (NX - 1) * 1.5, urcrnrlat=-60 + (NY - 1) * 1.5,
#                             resolution='i')
#        mapproj.drawcoastlines()
#        lonproj, latproj = mapproj(dx, dy)      #poject grid
#        CS1 = mapproj.contourf(lonproj, latproj, R_forec_residual_modes[k, :, :], clevs,
#                               cmap=barra, vmin=-1, vmax=1) #extended generate pretty colorbar
#        plt.title('Forec. resid. modes - ' + models[k].upper(), fontsize=22)
#    plt.suptitle('Correlation between Observations vs: ', fontsize=26, x=0.52, y=0.93)
#    plt.subplots_adjust(right=0.9, top=0.90)#, wspace=None, hspace=None)
#    cbar_ax = fig.add_axes([0.92, 0.33, 0.03, 0.35])
#    cax=fig.colorbar(CS1, cax=cbar_ax, orientation='vertical')
#    ticklabs = cax.ax.get_yticklabels()
#    cax.ax.set_yticklabels(ticklabs, fontsize=20)
#    plt.savefig(RUTA_OUT + 'skill_week' + str(i + 4) + '_bis.png', dpi=300,
#                bbox_inches='tight', orientation='portrait') 
#clevs = np.linspace(-2, 2, 11)
#barra = plt.cm.bwr #colorbar
#for i in [0]:
#    f = np.load('./outputs/skill_predic_' + WEEK[i] + '_bis.npz')
#    print(f.keys)
#    R_forec_predic_modes = f['arr_0']
#    print(np.min(R_forec_predic_modes))
#    print(np.max(R_forec_predic_modes))
#    print(R_forec_predic_modes)
#    R_forec_anom = f['arr_2']
#    dx = f['arr_4']
#    dy = f['arr_5']
#    fig = plt.figure(figsize=(30, 9), dpi=300)
#    for k in range(5):
#        plt.subplot(1, 5, k + 1)
#        mapproj = bm.Basemap(projection='cyl',llcrnrlon = 270, llcrnrlat = -60,
#                             urcrnrlon = 270+(NX-1)*1.5, urcrnrlat = -60+(NY-1)*1.5,
#                             resolution = 'i')
#        mapproj.drawcoastlines()
#        lonproj, latproj = mapproj(dx, dy)      #poject grid
#        CS1 = mapproj.contourf(lonproj, latproj, R_forec_predic_modes[k, :, :] /
#                               R_forec_anom[k, :, :], clevs,
#                               cmap=barra,vmin=-2,vmax=2, extend='both') #extended generate pretty colorbar
#        plt.title('Skill vs Predic modes - ' + models[k].upper(), fontsize=22)
#    plt.suptitle('Correlation between Observations vs: ', fontsize=26, x=0.52, y=0.93)
#    plt.subplots_adjust(right=0.9, top=0.90)#, wspace=None, hspace=None)
#    cbar_ax = fig.add_axes([0.92, 0.33, 0.028, 0.30])
#    cax=fig.colorbar(CS1, cax=cbar_ax, orientation='vertical')
#    ticklabs = cax.ax.get_yticklabels()
#    cax.ax.set_yticklabels(ticklabs, fontsize=20)
#    plt.savefig(RUTA_OUT + 'ratio_skill_predic_week' + str(i + 1) + '_bis.png', dpi=300,
#                bbox_inches='tight', orientation='portrait', papertype='A4') 
#for i in [0]:
#    f = np.load('./outputs/skill_predic_' + WEEK[i] + '_bis.npz')
#    R_forec_predic_modes = f['arr_0']
#    R_forec_residual_modes = f['arr_1']
#    R_forec_anom = f['arr_2']
#    R_obs_predic_modes = f['arr_3']
#    #    #f['R_obs_predic_modes1']
#    dx = f['arr_4']
#    dy = f['arr_5']
#    mapproj = bm.Basemap(projection='cyl',llcrnrlon = 270, llcrnrlat = -60,
#                         urcrnrlon = 270+(NX-1)*1.5, urcrnrlat = -60+(NY-1)*1.5,
#                         resolution = 'i')
#    lonproj, latproj = mapproj(dx, dy)      #poject grid
#    Array = np.zeros_like(dx, dtype=bool)
#    for j in range(dx.shape[0]):
#        for k in range(dx.shape[1]):
#            if mapproj.is_land(dx[j, k],dy[j, k]):
#                Array[j, k] = False
#            else:
#                Array[j, k] = True
#    aux = ma.masked_array(R_obs_predic_modes, mask=Array)
#    print('obs', np.reshape(aux, NX*NY).mean())
#    for k in range(5):
#        aux = ma.masked_array(R_forec_anom[k,:,:], mask=Array)
#        print(models[k], 'Anom', np.reshape(aux, NX*NY).mean())
#        aux = ma.masked_array(R_forec_predic_modes[k,:,:], mask=Array)
#        print(models[k], ' PM ', np.reshape(aux, NX*NY).mean())
#        aux = ma.masked_array(R_forec_residual_modes[k, :, :], mask=Array)
#        print(models[k], ' Res ', np.reshape(aux, NX*NY).mean())
#
fig = plt.figure(figsize=(30, 28), dpi=300)
colores = ['r', 'b', 'g', 'y', 'k']
for i in range(4):
    f = np.load('./outputs/pss_exp_var_' + WEEK[i] + '.npz')
    pss = f['arr_2']
    obs_var = f['arr_1']
#    print(obs_var*100)
    plt.subplot(2, 2, i + 1)
    for k in range(5):
        print(models[k], pss[k, 0:])

        CS = plt.scatter(pss[k, :], obs_var[0:7] * 100, color=colores[k], label=models[k].upper(), s=360)
    plt.xlabel('PSS', fontsize=24)
    plt.tick_params(axis='both', which='major', labelsize=18)
    plt.xlim([0, 1])
    plt.ylim([1,45])
    plt.ylabel('Observed explained variance', fontsize=24)
    plt.title('Week ' + str(i + 1), fontsize=28)
plt.subplots_adjust(right=0.9, top=0.90)
plt.legend(bbox_to_anchor=(1.05, 0), loc='lower left', borderaxespad=0., fontsize=24)
#plt.savefig(RUTA_OUT + 'pss_obs_var_new.png', dpi=300,  bbox_inches='tight')
##
