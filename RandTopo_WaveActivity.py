import numpy as np
import os
from scipy.io import netcdf as scipy_netcdf
from matplotlib import pyplot as plt
from scipy.signal import firwin, filtfilt, lfilter
from scipy.fftpack import fft, fftfreq, fftshift, ifft, fft2, ifft2
import MontgomeryPot as MP

def computeDivergence(vectorX,vectorY, deltaX, deltaY, layer):
	
	gradientX = np.gradient(vectorX[layer,:,:], deltaY,deltaX)
	gradientY = np.gradient(vectorY[layer,:,:], deltaY,deltaX)
	
	
	#return gradientX[1] + gradientY[0]
	return gradientX[1]
meanDepth = 4000.0
G0   = 0.980 
RHO0 = 1025.0 
reducedGravity = np.asarray([0.005,0.005,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025,0.0025])
densityLayers = np.zeros(reducedGravity.shape[0],dtype='float')
refThickness  = np.asarray([150.0,150.0,275.0,275.0,275.0,275.0,650.0,650.0,650.0,650.0])
for iLayer in range(0,reducedGravity.shape[0]):
	if iLayer == 0:
		densityLayers[iLayer] = RHO0
	else:
		densityLayers[iLayer] = (reducedGravity[iLayer-1]/9.80) * RHO0 + densityLayers[iLayer-1]



topoName = 'topo'

iDir = 5
dirName = '/home/chris/GOLD/RandTopo_Smooth_Processing/'
zonalWaveActivityFileName = dirName + 'waveActivity_transient_RandTopoSmooth.nc'
print zonalWaveActivityFileName	



topoFileName  = dirName + 'isoTopo_' + str(iDir) + '.nc'
EKE_fileName = dirName + 'postprocessed_IsoTopo' + str(iDir) + '.nc'
print EKE_fileName

waveFile = scipy_netcdf.netcdf_file(zonalWaveActivityFileName, 'r')
topoFile  = scipy_netcdf.netcdf_file(topoFileName, 'r')
ekeFile  = scipy_netcdf.netcdf_file(EKE_fileName, 'r')

topo5        = meanDepth - topoFile.variables[topoName][:,:]

#Get dimensions
xName = 'x'
yName = 'y'
zlName = 'zl'
ekeName  = 'eke'


x  = waveFile.variables[xName][:]
y  = waveFile.variables[yName][:]
zl = waveFile.variables[zlName][:]
topoFile.close()

deltaX = x[1]-x[0]
deltaY = y[1]-y[0]


EPVector_5 = []


EPVector_5.append(waveFile.variables['EP_X'][:,:,:])
EPVector_5.append(waveFile.variables['EP_Y'][:,:,:])
EPVector_5.append(waveFile.variables['EP_Z'][:,:,:])

ekeTimeMean5 = ekeFile.variables[ekeName][:,:,:]
vTimeMean5 = ekeFile.variables['v_time_mean'][:,:,:]
layerToPlot = 4


print vTimeMean5.shape
print ekeTimeMean5.shape
print x.shape
print y.shape
ekeFile.close()
waveFile.close()
