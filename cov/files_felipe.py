# # Austral summer sub-seasonal precipitation predictable modes over South America #
# 
# ## Questions to be addressed: ##
# + How do the analyzed S2S models (and the multi-model ensemble of these models) represent the dominant modes of sub-seasonal precipitation variability during the extended austral summer at different lead times? 
# - What modes of sub-seasonal precipitation variability contribute to the potential predictability during the austral summer? 
# * What percentage of the actual S2S models prediction ability comes from the identified modes of sub-seasonal precipitation variability? 
# 

# In[2]:
import numpy as np
import xarray as xr
import pandas as pd
from scipy import stats
covariance = True
WEEK = ['w1', 'w2', 'w3', 'w4']
RUTA = './outputs/'
RUTA_OUT = './figures_caio/'
NWEEKS = 242
NX = 41
NY = 51
NMODES = 4
NPOINTS = NX*NY
models = ['cma', 'eccc', 'ecmwf', 'ncep', 'mmem']
#np.savez('./outputs/EOF_analsys_w1.npz',obs_pc_w1=obs_pc,obs_evec_w1=obs_evec,
#         models_pc_w1=models_pc_sorted, models_evec_w1=models_evec_sorted, dx=dx, dy=dy)
tiempo = []
for i in range(1999,2011):
    time = pd.date_range(str(i)  + '-01-05', freq='7D', periods= 52).tolist()
    tiempo.append(time)

tiempo = [y for x in tiempo for y in x]

tiempo = pd.DatetimeIndex(tiempo)

ds = xr.DataArray(np.arange(len(tiempo)), coords=[('time', tiempo)])

tiempos = ds.time.sel(**{'time':slice('1999-11-01', '2010-03-31')})

tiempos = tiempos['time'].sel(time=np.logical_or(tiempos.time['time.month']>=11, tiempos.time['time.month']<=3))

mode = np.arange(1,5)

PCs = np.empty([4, 5, 4, 242])
PCs_obs = np.empty([4, 4, 242])

for i in range(4):
    FILE = RUTA + 'EOF_analsys_' + WEEK[i] + '.npz'
    f = np.load(FILE)
    for j in range(4):
        signo = -1 if np.logical_and(j == 2, i==0) else 1

        PCs_obs[i, j, :] = signo * f['obs_pc_' + WEEK[i]][j, :]

        for k in range(5):
            PCs[i, k, j, :] = signo * f['models_pc_' + WEEK[i]][k, j, :]


ds = xr.DataArray(PCs, coords=[('week', WEEK), ('model', models),('mode', mode), ('IC', tiempos)])
ds.name = 'PC'
ds.to_netcdf('models_pc.nc4')
ds1 = xr.DataArray(PCs_obs, coords=[('week', WEEK), ('mode', mode), ('IC', tiempos)])
ds1.name = 'PC'
ds1.to_netcdf('obs_pc.nc4')


#figpd.DatetimeIndex(tiempo)
#for i in [3]:
##    fig = plt.figure(figsize=(26, 42), dpi=300)
#    for j in range(4):

#
#for i in [3]:#range(3):
##    fig = plt.figure(figsize=(41, 10), dpi=300)  #fig size in inches
#    for j in np.arange(4):
#        signo = -1 if np.logical_and(j == 2, i==0) else 1 
#        FILE = RUTA + 'EOF_analsys_' + WEEK[j] + '.npz'
#        f = np.load(FILE)
#        y_obs = signo * f['obs_pc_' + WEEK[j]][i, :]
##        ax.plot(x, y, color=c, linewidth=1.5, label='Obs')
#        for k in np.arange(5):
##            c = next(color)
#            y = #            print(models[k], str(j+1), np.corrcoef(y,y_obs)[0,1])
#
#
#
