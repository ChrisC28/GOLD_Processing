import numpy as np
import os
from scipy.io import netcdf as scipy_netcdf
from matplotlib import pyplot as plt
from scipy.signal import firwin, filtfilt, lfilter
from scipy.fftpack import fft, fftfreq, fftshift, ifft, fft2, ifft2
import PeakDetection as PD


meanDepth = 4000.0
iLayer  = 1

outputDirectory = '/home/chris/GOLD/IsoTopo_1_Processing/'

modelOutputFileName = 'postprocessed_IsoTopo1.nc'
topoFileName   = 'isoTopo1.nc'

modelFile = scipy_netcdf.netcdf_file(outputDirectory + modelOutputFileName, 'r')
topoFile  = scipy_netcdf.netcdf_file(outputDirectory + topoFileName, 'r')
#Get dimensions
xName = 'x'
yName = 'y'
zlName = 'zl'
ziName = 'zi'

x  = modelFile.variables[xName][:]
y  = modelFile.variables[yName][:]
zl = modelFile.variables[zlName][:]
zi = modelFile.variables[ziName][:]

deltaX = x[1]-x[0]
deltaY = y[1]-y[0]


vName    = 'v_time_mean'
uName    = 'u_time_mean'

ekeName  = 'eke'
topoName = 'topo'

uTimeMean   = modelFile.variables[uName][:,:,:]
vTimeMean   = modelFile.variables[vName][:,:,:]
ekeTimeMean = modelFile.variables[ekeName][:,:,:]
topo        = meanDepth - topoFile.variables[topoName][:,:]

modelFile.close()
topoFile.close()



