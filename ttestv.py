#!/usr/bin/env python3
############################################
######## author: Georgia Pantelidou ########
######### created on January 2021 ##########
############################################

import netCDF4
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from scipy import stats

#import data
data = netCDF4.Dataset(r'C:\Users\gpant\OneDrive\Desktop\scripts\MIROC\1800-1849_tas_mean_lonlat.nc')
data_base = netCDF4.Dataset(r'C:\Users\gpant\OneDrive\Desktop\scripts\MIROC\950-999_tas_base_mean_lonlat.nc')
data_anom = netCDF4.Dataset(r'C:\Users\gpant\OneDrive\Desktop\scripts\MIROC\1800-1849_anom_lonlat.nc')
#create variables
lat = data.variables['lat'][:] 
lon = data.variables['lon'][:]                   
time = data.variables['time'][:]
tas = data.variables['tas'][:]
tas_base = data_base.variables['tas'][:]
tas_anom = data_anom.variables['tas'][:]



n = len(time) #sample size
#with axis=0 we run in the time dimension and calculate temprerature means and variances 
mean_tas = tas.mean(axis=0)
mean_tas_base = tas_base.mean(axis=0)
var_tas = tas.var(axis=0)
var_tas_base = tas_base.var(axis=0)
#two-sample independent t-test
tvalue = np.abs(mean_tas - mean_tas_base)/np.sqrt((var_tas + var_tas_base)/(n-1))
df = 2*n - 2 #degrees of freedom
pvalue = (1 - stats.t.cdf(tvalue,df=df))*2


#res_list = [] 
#for i in range(0, len(lon)) : 
 #   if lon[i] > 20 and lon[i] < 25 : 
  #      res_list.append(i) 
#print(res_list)



print(tvalue)
print("\n")
print(pvalue)
print("\n")

lat_max = np.max(lat)
lat_min =np.min(lat)
lon_min = np.min(lon)
lon_max = np.max(lon)

mp = Basemap(projection= 'mill', 
             llcrnrlon = lon_min,
             llcrnrlat = lat_min,
             urcrnrlon = lon_max,
             urcrnrlat = lat_max,
             resolution = 'l')

lon,lat = np.meshgrid(lon,lat)
x,y = mp(lon,lat)

#if we receive a p-value < 0.05, we can reject the null hypothesis and state that there is a significant difference.
#stat_sign = np.zeros((14,25))
    
#for lat in range(14):
 #   for lon in range(25):
  #      if pvalue[lat, lon] < 0.05:
   #         stat_sign[lat, lon] = pvalue[lat, lon]  
#print(stat_sign)

#contour levels
colorbar_levels = np.linspace(start=-1, stop=1.0, num=17) #17,9,10
#colorbar_levels = np.linspace(start=-1.5, stop=1.5, num=27) #17,9,10
#colorbar_levels = np.linspace(start=-1.5, stop=1.5, num=13) #17,9
#create a map of Europe with temperarure anomalies
europe_map = mp.contourf(x,y, np.squeeze(tas_anom), cmap = 'RdBu_r', levels=colorbar_levels)

#create a second map which overlays the map of europe with statistical significant points


#levels = [stat_sign[:,:].min(), 0.000000000001, stat_sign[:,:].max()]
#mp.contourf(x, y, np.squeeze(stat_sign)[:,:], cmap = 'RdBu_r', alpha=0, hatches=["","."])

levels = [pvalue[:,:].min(), 0.05, pvalue[:,:].max()]
mp.contourf(x, y, pvalue[:,:], levels=levels, alpha=0, hatches=[".",""])



mp.drawcoastlines(zorder=1)
#mp.drawcountries()
mp.drawmapboundary(fill_color='white')
cbar = mp.colorbar(europe_map, location = 'right' , pad = '5%')
cbar.ax.tick_params(labelsize=10) 


plt.title('MIROC-ES2L  1800-1849 (950-999)', fontsize=18)
cbar.minorticks_on()
plt.show()
