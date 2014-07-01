import numpy as np
import os
from scipy.io import netcdf as scipy_netcdf
from matplotlib import pyplot as plt
from scipy.signal import firwin, filtfilt, lfilter
from scipy.fftpack import fft, fftfreq, fftshift, ifft, fft2, ifft2

def FilterField2D(inputField,cutoffFreqX,cutoffFreqY,samplingIntervalX,samplingIntervalY):
	
	[nY,nX] = inputField.shape
	
	kernalLength = 20.0
	
	nyquistFrequencyX = 0.5/samplingIntervalX
	nyquistFrequencyY = 0.5/samplingIntervalY
	
	tapWeightsX = firwin(kernalLength,cutoffFreqX,window='blackman',nyq=nyquistFrequencyX)
	tapWeightsY = firwin(kernalLength,cutoffFreqY,window='blackman',nyq=nyquistFrequencyY)

	filteredField = np.zeros([nY,nX],dtype='float64')
	for iY in range(0,nY):
		filteredField[iY,:] = filtfilt(tapWeightsX,[1.0],inputField[iY,:])
	
	for iX in range(0,nX):
		filteredField[:,iX] = filtfilt(tapWeightsY,[1.0],inputField[:,iX])

	return filteredField
#Get Data

modelDataFileName = '/home/chris/GOLD/postprocessedOut.nc'
topoFileName = '/home/chris/ETOPO1/randTopo_01.nc'


dataFile  = scipy_netcdf.netcdf_file(modelDataFileName,'r')
topoFile  = scipy_netcdf.netcdf_file(topoFileName,'r')

#Get variables


x  = dataFile.variables['x'][:]
y  = dataFile.variables['y'][:]
zl = dataFile.variables['zl'][:]
zi = dataFile.variables['zi'][:]

uTimeMean = dataFile.variables['u_time_mean'][:,:,:]
vTimeMean = dataFile.variables['v_time_mean'][:,:,:]
ekeTimeMean = dataFile.variables['eke'][:,:,:]

topo = topoFile.variables['topo'][:,:]


dataFile.close()
topoFile.close()

cutoffPeriodX = 100.0 * (x[1]-x[0])
cutoffPeriodY = 50.0 * (y[1]-y[0])
LPFilteredTopo = FilterField2D(topo,1.0/cutoffPeriodX,1.0/cutoffPeriodY,x[1]-x[0],y[1]-y[0])



fig = plt.figure(1)
ax  = fig.add_subplot(1,1,1)
ax.contourf(x,y,uTimeMean[0,:,:],15,cmap=plt.cm.jet)
ax.contour(x,y,LPFilteredTopo[1::,:],10,colors='k')
plt.show()
