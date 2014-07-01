import numpy as np
import os
from scipy.io import netcdf as scipy_netcdf
from matplotlib import pyplot as plt
from scipy.signal import firwin, filtfilt, lfilter
from scipy.fftpack import fft, fftfreq, fftshift, ifft, fft2, ifft2
import PeakDetection as PD

def FindPeaksSpectrum(inputFFT):
	
	nPoints = inputFFT.shape[-1]
	inputSpectrum = inputFFT * inputFFT.conj()
	firstDeriv = np.diff(inputSpectrum,axis=-1)
	zero_crossings = np.where(np.diff(np.sign(firstDeriv[0,nPoints/2:])))[0]
	
	
	for iCrossing in range(len(zero_crossings)):
		zero_crossings[iCrossing] = zero_crossings[iCrossing] + 1 + nPoints/2
	
	
	return zero_crossings
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

meanDepth = 4000.0
iLayer  = 2

outputDirectory = '/home/chris/GOLD/IsoTopo_4_Processing/'

modelOutputFileName = 'postprocessed_IsoTopo4.nc'
topoFileName   = 'isoTopo_4.nc'

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

#Average over a latitude band
#lowerLatitude  = 500.0
#upperLatitude = 1250.0

lowerLatitude  = 700.0
upperLatitude = 1350.0


lowerIndex = np.nonzero(y >= lowerLatitude)[0][0]
upperIndex = np.nonzero(y >= upperLatitude)[0][0]
averagedEKE = np.mean(ekeTimeMean[:,lowerIndex:upperIndex,:],axis=1)
averagedV   = np.mean(vTimeMean[:,lowerIndex:upperIndex,:],axis=1)

peakEKE, maxIndex, endIndex, integratedEKE = PD.StormTrackMetrics(averagedEKE[iLayer,:])

if maxIndex < endIndex:
	stormTrackLength = (endIndex - maxIndex) * deltaX
else:
	stormTrackLength = ((x.shape[0]-maxIndex) + (endIndex) ) * deltaX 


print peakEKE
print "storm track length:", stormTrackLength
print "integrated EKE:", integratedEKE
print x[maxIndex]
print x[endIndex]
print x[endIndex] - x[maxIndex]

averagedV_FT = fftshift(fft(averagedV[:,maxIndex::],axis=1))
x_freq      = fftshift(fftfreq(x[maxIndex::].shape[0],d=deltaX))

print "average V FT:", averagedV_FT.shape
fftPeakIndex = FindPeaksSpectrum(averagedV_FT)

peakWavenumber  = x_freq[fftPeakIndex]

peakWavelength  = 1.0/peakWavenumber
print "Peak Wavelength:", peakWavelength
print "PeakWaveNumber:", peakWavenumber

print "PeakWaveNumber:", peakWavenumber*deltaX*x.shape[0]



fig = plt.figure(0)
ax  = fig.add_subplot(1,1,1)
ax.set_title("Time Mean Zonal Velocity (m/s)")
ax.contourf(x,y,uTimeMean[iLayer,:,:],15,cmap=plt.cm.jet)
ax.contour(x,y,topo[1::,:],10,colors='k')
ax.plot(x,y[upperIndex]*np.ones(x.shape[0],dtype='float64'),'k')
ax.plot(x,y[lowerIndex]*np.ones(x.shape[0],dtype='float64'),'k')
ax.set_xlabel('x (km)')
ax.set_ylabel('y (km)')


fig = plt.figure(1)
ax  = fig.add_subplot(1,1,1)
ax.set_title("Time Mean Eddy Kinetic Energy (m2/s2)")
ax.contourf(x,y,ekeTimeMean[iLayer,:,:],15,cmap=plt.cm.jet)
ax.contour(x,y,topo[1::,:],10,colors='k')
ax.plot(x,y[upperIndex]*np.ones(x.shape[0],dtype='float64'),'k')
ax.plot(x,y[lowerIndex]*np.ones(x.shape[0],dtype='float64'),'k')
ax.set_xlabel('x (km)')
ax.set_ylabel('y (km)')





fig = plt.figure(2)
ax  = fig.add_subplot(1,1,1)
ax.set_title("Time Mean Meridional Velocity (m/s)")
ax.contourf(x,y,vTimeMean[iLayer,:,:],15,cmap=plt.cm.jet)
ax.contour(x,y,topo[1::,:],10,colors='k')
ax.plot(x,y[upperIndex]*np.ones(x.shape[0],dtype='float64'),'k')
ax.plot(x,y[lowerIndex]*np.ones(x.shape[0],dtype='float64'),'k')
ax.set_ylabel('y (km)')
ax.set_xlabel('x (km)')



fig = plt.figure(3)
ax  = fig.add_subplot(2,1,2)
ln1=ax.plot(x,averagedEKE[iLayer,:]*1000.0,'k',label='EKE')
ax.plot([x[maxIndex],x[maxIndex]],[np.min(averagedEKE[iLayer,:])-0.01,40.0],'k:')
ln3=ax.plot([x[endIndex],x[endIndex]],[np.min(averagedEKE[iLayer,:])-0.01,40.0],'k:')
ax2 = ax.twinx()
ln2=ax2.plot(x,topo[(upperIndex),:],'k--',label='Topo')
ax2.set_ylabel('Topo (m)')
lns = ln1+ln2
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc='upper left', shadow=True, prop={'size':10})


ax.set_xlim([0,x[-1]])
ax.set_ylim([0,np.max(averagedEKE[iLayer,:])])
ax.set_yticks(np.arange(0,(0.041)*1000.0,(0.01)*1000.0))
ax.annotate('b',xy=(0.01,1.025),xycoords='axes fraction',fontsize=15)
ax.set_xlabel('x (km)')

ax.set_ylabel('EKE (J.m$^{-3}$)')


ax  = fig.add_subplot(2,1,1)
ln1=ax.plot(x,averagedV[iLayer,:],'k',label='V')
ax2 = ax.twinx()
ln2=ax2.plot(x,topo[(upperIndex),:],'k--',label='Topo')
ax2.set_ylabel('Topo (m)')
lns = ln1+ln2
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc='upper left', shadow=True, prop={'size':10})
ax.annotate('a',xy=(0.01,1.025),xycoords='axes fraction',fontsize=15)
ax.set_xlim([0,x[-1]])
ax.set_ylabel('V (m.s$^{-1}$)')
ax.get_xaxis().set_visible(False)


'''
ax  = fig.add_subplot(3,1,3)
ax.plot(x,topo[(upperIndex),:],'k')
ax.set_xlabel('x (km)')
ax.set_ylabel('Topo (m)')
ax.annotate('c',xy=(0.01,1.025),xycoords='axes fraction',fontsize=15)
ax.set_xlim([0,x[-1]])

fig = plt.figure(4)
ax  = fig.add_subplot(1,1,1)
ax.plot(x_freq[x_freq.shape[0]/2::],(averagedV_FT*averagedV_FT.conj())[1,x_freq.shape[0]/2::])
'''


iLayer = 0
fig = plt.figure(5)
ax  = fig.add_subplot(3,1,1)
#ax.set_title("Time Mean Zonal Velocity (m/s)")
uVals = np.arange(-0.1,0.51,0.025)
uValsCbar = np.arange(-0.1,0.505,0.1)
thicknessAveU = np.zeros((uTimeMean[1,...]).shape,dtype='float64')
for iLayer in range(0,zl.shape[0]):
	thicknessAveU[:,:] = thicknessAveU[:,:] + (refThickness[2.0*iLayer] * uTimeMean[iLayer,:,:])
thicknessAveU = thicknessAveU / refThickness[::2].sum()


cs = ax.contourf(x,y,thicknessAveU,uVals,cmap=plt.cm.RdBu_r)
ax.annotate('a',xy=(0.01,1.025),xycoords='axes fraction',fontsize=15)


cbar = fig.colorbar(cs,ax=ax, ticks=uValsCbar)
cbar.ax.set_title('m/s')

ax.contour(x,y,topo[1::,:],10,colors='k')
#ax.plot(x,y[upperIndex]*np.ones(x.shape[0],dtype='float64'),'k')
#ax.plot(x,y[lowerIndex]*np.ones(x.shape[0],dtype='float64'),'k')
ax.set_ylabel('y (km)')
ax.get_xaxis().set_visible(False)




thicknessAveV = np.zeros((vTimeMean[1,...]).shape,dtype='float64')
for iLayer in range(0,zl.shape[0]):
	thicknessAveV[:,:] = thicknessAveV[:,:] + (refThickness[2.0*iLayer] * vTimeMean[iLayer,:,:])
thicknessAveV = thicknessAveV / refThickness[2::].sum()


ax  = fig.add_subplot(3,1,2)
vVals = np.arange(-0.15,0.151,0.01)
vValsCbar = np.arange(-0.15,0.151,0.05)

cs = ax.contourf(x,y,thicknessAveV[:,:],vVals,cmap=plt.cm.RdBu_r)
cbar = fig.colorbar(cs,ax=ax,ticks = vValsCbar)
cbar.ax.set_title('m/s')
#cbar.ax.set_yticklabels(['{0:.0e}'.format( -0.08),'{0:.0e}'.format( -0.04),0,'{0:.0e}'.format( 0.04),'{0:.0e}'.format( 0.08) ]) 
ax.annotate('b',xy=(0.01,1.025),xycoords='axes fraction',fontsize=15)

ax.contour(x,y,topo[1::,:],10,colors='k')

#ax.plot(x,y[upperIndex]*np.ones(x.shape[0],dtype='float64'),'k')
#ax.plot(x,y[lowerIndex]*np.ones(x.shape[0],dtype='float64'),'k')
ax.set_ylabel('y (km)')
ax.get_xaxis().set_visible(False)


thicknessAveEKE = np.zeros((ekeTimeMean[1,...]).shape,dtype='float64')
for iLayer in range(0,zl.shape[0]):
	thicknessAveEKE[:,:] = thicknessAveEKE[:,:] + (densityLayers[2.0*iLayer]*refThickness[2.0*iLayer] * ekeTimeMean[iLayer,:,:])
thicknessAveEKE = thicknessAveEKE / refThickness[::2].sum()

ax  = fig.add_subplot(3,1,3)
#ekeVals = np.arange(0,0.20+0.01,0.025) *densityLayers[iLayer/2]
ekeVals = np.arange(0,51,5)
#ekeValsCbar = np.arange(0,0.0+0.01,0.05)* densityLayers[iLayer/2]
ekeValsCbar = [0.0,10,20.0,30.0,40.0,50]
cs = ax.contourf(x,y,thicknessAveEKE,ekeVals,cmap=plt.cm.RdBu_r)
ax.annotate('c',xy=(0.01,1.025),xycoords='axes fraction',fontsize=15)


cbar = fig.colorbar(cs,ax=ax,ticks=ekeValsCbar)
cbar.ax.set_title('J.m$^{-3}$')

ax.contour(x,y,topo[1::,:],10,colors='k')

#ax.plot(x,y[upperIndex]*np.ones(x.shape[0],dtype='float64'),'k')
#ax.plot(x,y[lowerIndex]*np.ones(x.shape[0],dtype='float64'),'k')
ax.set_xlabel('x (km)')
ax.set_ylabel('y (km)')


xIndex3000 = np.nonzero(x >= 2000.0)[0][0]


fig = plt.figure(6)
ax = fig.add_subplot(1,1,1)
ax.plot(y,uTimeMean[iLayer,:,xIndex3000],'k')


plt.show()
















