import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from netCDF4 import Dataset
from scipy.interpolate import griddata
from cmocean import cm
from matplotlib.colors import LightSource
import gsw
from matplotlib import ticker
from datetime import datetime

# Define the variables
vars = ['T', 'C', 'S', 'chla', 'CDOM', 'BB', 'DO']
varN = ['Temperature', 'Conductivity', 'Salinity', 'Chlorophyll a', 'CDOM', 'Backscatter', 'Dissolved Oxygen']
units = ['^oC', 'S m^-1', 'PSU']
clims = [[21, 29], [4, 7], [39, 40.5], [0, 1], [0.2, 0.7], [0, 3e-4], [20, 180]]

cbars = [cm.thermal, cm.solar, cm.haline, cm.algae, cm.algae, cm.matter, cm.oxy]


# Read the netCDF file
ncfile = Dataset('slocum_timeseries.nc')

# Read the variables
Time = ncfile.variables['time'][:] / (3600 * 24) + 719529  # 719529 is the difference in days between 1970 and 0000 in proleptic Gregorian calendar
T = ncfile.variables['temperature'][:]
C = ncfile.variables['conductivity'][:]
S = ncfile.variables['salinity'][:]
P = ncfile.variables['press'][:]
depth = ncfile.variables['depth'][:]
chla = ncfile.variables['chla'][:]
DO = ncfile.variables['DO'][:]
CDOM = ncfile.variables['CDOM'][:]
BB = ncfile.variables['BB'][:]
Dens = ncfile.variables['dens'][:]
Lat = ncfile.variables['lat'][:]
Lon = ncfile.variables['lon'][:]

# Convert time to datetime
dates = mdates.num2date(Time)

# Calculate time step
Timestep = (np.max(Time) - np.min(Time)) / 1000

# Create time bins
TimeBins = np.arange(np.min(Time), np.max(Time), Timestep)

# Define depths
Depths = np.arange(-1, -601, -1)

# Generate a grid for interpolation
XX, YY = np.meshgrid(TimeBins, Depths)

# Find indices to remove
ind = np.where((depth > 0) | np.isnan(depth) | (depth < -1000))

# Remove indices from data arrays
depth = np.delete(depth, ind)
Time = np.delete(Time, ind)
T = np.delete(T, ind)
C = np.delete(C, ind)
S = np.delete(S, ind)
P = np.delete(P, ind)
chla = np.delete(chla, ind)
DO = np.delete(DO, ind)
CDOM = np.delete(CDOM, ind)
BB = np.delete(BB, ind)
Dens = np.delete(Dens, ind)
Lat = np.delete(Lat, ind)
Lon = np.delete(Lon, ind)

# Interpolate to get data on regular grid
Dbin = griddata((Time, depth), Dens, (XX, YY), method='linear')


# Define the custom formatting function
def format_date(value, tick_number):
    date_value = datetime.fromtimestamp(value)
    formatted_date = date_value.strftime('%Y:%m:%d %H:%M:%S')
    return formatted_date

def plotHovm(lm, tm, Vm, Dens, climits, cLabel, Title, cbar, filename):
    DensLevels = [1025, 1026, 1027, 1028, 1028.5]
    fig, ax = plt.subplots(figsize=(10, 7.5))
    
    pcolor = ax.pcolormesh(lm, tm, Vm, cmap=cbar, shading='auto')
    if not np.isnan(climits).any():
        pcolor.set_clim(climits)
    cb = plt.colorbar(pcolor)
    cb.set_label(cLabel)
    
    contour = ax.contour(lm, tm, Dens, DensLevels, linewidths=2, colors='w')
    ax.clabel(contour, inline=True, fontsize=10, colors='w')

    ax.set_ylim([-600, 0])
    ax.set_title(Title)
    
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(format_date))
    fig.autofmt_xdate()
    
    plt.savefig(filename)

def TSplot(T, S, depth, filename):
    fig, ax = plt.subplots(figsize=(10, 7.5))
    plt.scatter(S, T, c=-depth, marker='.', cmap=cm.balance)
    plt.colorbar(label="Depth")
    plt.clim(0, 600)
    
    SS = np.arange(38.6, 41.001, 0.001)
    TT = np.arange(21, 29.01, 0.01)
    SP, TP = np.meshgrid(SS, TT)

    PP = np.zeros(SP.shape)
    D = gsw.density.rho(SP, TP, PP)
    
    plt.contour(SP, TP, D, colors='#666666')
    plt.xlabel("Salinity (PSU)")
    plt.ylabel("Temperature (^oC)")
    plt.xlim([38.5, 41])
    
    plt.savefig(filename)



vars  = ['T', 'C', 'S', 'chla', 'CDOM', 'BB', 'DO']
varN  = ['Temperature', 'conductivity', 'Salinity', 'chlorophyl a', 'CDOM', 'Backscatter', 'Dissolved Oxygen']
# cbars = ['thermal', 'thermal', 'haline', 'algae', 'matter', 'matter', 'parula']
outfolder = './png/'

for i in range(len(vars)):
    var = vars[i]
    F = griddata((Time, depth), globals()[var], (XX, YY), method='linear', fill_value=np.nan)
    plotHovm(XX, YY, F, Dbin, clims[i], varN[i], varN[i], cbars[i], outfolder + varN[i] + '.png')

TSplot(T, S, depth, outfolder + 'TS' + '.png')
