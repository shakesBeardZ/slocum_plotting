
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks
from numpy import linspace, meshgrid
from scipy.interpolate import griddata
from typing import Tuple
import cmocean 

def create_mask(X, Y, timeData, zData, peaks):
    mask = np.ma.make_mask(None)
    maxDepth = np.min(-zData)
    yPeaks = -zData[peaks]
    xPeaks = timeData[peaks]

    length = len(xPeaks)
    for i in range(length):
        mask_condition = (Y > maxDepth) & (Y < yPeaks[i]) & (X >= xPeaks[i]) & (X < xPeaks[i] + 1)
        mask = np.logical_or(mask, mask_condition)

    return mask

def plot_contourf(X, Y, Z, cmap):
    return plt.contourf(X, Y, Z, cmap=cmap)

def plot_colorbar(img, dataName, dataUnit):
    cbar = plt.colorbar(img)
    cbar.ax.get_yaxis().labelpad = 15
    cbar.ax.set_ylabel(f'{dataName} [{dataUnit}]', rotation=270)
    return cbar

def plotImageTimeSeries(data, dataUnit, dataName, zData, zUnit, timeData,
              cmap, filename, targetDir):
    filename = "_".join(filename.split("_")[:-1])
    fname = f'{targetDir}/{filename}_timeline.png'
    print (f'plotImage {filename} : len(data): {len(data)}')

    # Calculate plot limits
    dmin = np.min(-zData)
    xmax = np.nanmax(timeData)
    xmin = np.nanmin(timeData)
    onePerc = (xmax - xmin) / 100

    # Create figure and axes
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1)

    # Hide axes spines
    for spine in ['top', 'bottom', 'right', 'left']:
        ax.spines[spine].set_visible(False)

    # Set up axes
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.ylim(dmin-10, 0)
    plt.xlim(xmin-onePerc, xmax+onePerc)
    plt.yticks(fontsize=14)
    plt.xticks(fontsize=14)
    plt.tick_params(axis="both", which="both", bottom="off", top="off",
                    labelbottom="on", left="off", right="off", labelleft="on")

    # Create grid
    X, Y, Z = create_grid(timeData, -zData, data)

    # Create mask to handle incorrect interpolation
    maxDepth = np.min(-zData)
    peaks, _ = find_peaks(-zData, height=(maxDepth, -30), width=2)
    mask = create_mask(X, Y, timeData, zData, peaks)
    Z = np.ma.masked_array(Z, mask=mask)

    # Plot contourf, colorbar, and save the plot
    img = plot_contourf(X, Y, Z, cmap)
    plot_colorbar(img, dataName, dataUnit)
    ax.grid()
    plt.ylabel(f'Depth [{zUnit}]', fontsize=14)
    plt.xlabel('Time', fontsize=14)
    plt.title(dataName, fontsize=14, ha="center")
    plt.savefig(fname)
    plt.close()



def create_grid(x: np.ndarray, y: np.ndarray, z: np.ndarray, resX: int = 100, resY: int = 100) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Convert 3 column data to matplotlib grid.

    Parameters:
    x (np.ndarray): 1-D array representing the x-coordinates of the data points.
    y (np.ndarray): 1-D array representing the y-coordinates of the data points.
    z (np.ndarray): 1-D array representing the values at each data point.
    resX (int): Desired resolution in the x-direction. Default is 100.
    resY (int): Desired resolution in the y-direction. Default is 100.

    Returns:
    Tuple[np.ndarray, np.ndarray, np.ndarray]: Three 2-D arrays representing the x-coordinates, y-coordinates, and values at each point in the output grid.
    """
    xi = linspace(min(x), max(x), resX)
    yi = linspace(min(y), max(y), resY)
    X, Y = meshgrid(xi, yi)
    Z = griddata((x, y),z,(xi[None,:], yi[:,None]))
    
    return X, Y, Z


def readNCFile(filename):
    dataset = Dataset(filename, mode='r')

    variables = {
        'time': 'time',
        'temperature': 'temperature',
        'latitude': 'lat',
        'longitude': 'lon',
        'dissolved_oxygen': 'DO',
        'ctd_data_point_dive_number': '',
        'depth': 'depth',
        'nav_depth': 'navDepth',
        'salinity': 'salinity',
        'density': 'dens',
        'backscatter': 'BB',
        'CDOM': 'CDOM',
        'chlorophyl_a': 'chla',
        'conductivity': 'conductivity',
        'pressure': 'press'
    }

    data = {}
    for var, attribute in variables.items():
        try:
            data[var] = dataset.variables[attribute][:]
            data[f'{var}_units'] = dataset.variables[attribute].units
        except KeyError:
            print(f"Warning: {var} not found in the dataset.")
            data[var] = None
            data[f'{var}_units'] = None

    # Special handling for density
    # if data['density'] is not None:
    #     data['density'] -= 1000

    dataset.close()

    return data

def main():
    print(f"[INFO]: Plotting Slocume Glider data")
    data = readNCFile('unit_907-2023-133-4_timeseries.nc')
    plotImageTimeSeries(data['temperature'], data['temperature_units'], 'Temperature', data['nav_depth'], data['nav_depth_units'], data['time'],
        cmocean.cm.thermal, 'temperature', 'target')

if __name__ == '__main__':
    main()