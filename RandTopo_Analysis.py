import numpy as np
import os
from scipy.io import netcdf as scipy_netcdf
from matplotlib import pyplot as plt
from scipy.signal import firwin, filtfilt, lfilter
from scipy.fftpack import fft, fftfreq, fftshift, ifft, fft2, ifft2
import PeakDetection as PD

F0 = -1.0e-4
BETA = 1.5e-11
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


def FindPeaksSpectrum(inputFFT):
	
	nPoints = inputFFT.shape[-1]
	inputSpectrum = inputFFT * inputFFT.conj()
	firstDeriv = np.diff(inputSpectrum,axis=-1)
	zero_crossings = np.where(np.diff(np.sign(firstDeriv[0,nPoints/2:])))[0]
	
	
	for iCrossing in range(len(zero_crossings)):
		zero_crossings[iCrossing] = zero_crossings[iCrossing] + 1 + nPoints/2
	
	
	return zero_crossings

meanDepth = 4000.0

outputDirectory = '/home/chris/GOLD/RandTopo_Smooth_Processing/'

modelOutputFileName = 'postprocessed_RandSmooth.nc'
topoFileName   = 'randTopo_01_smooth.nc'

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


iLayer = 2
fig = plt.figure(0)
ax  = fig.add_subplot(3,1,1)
#ax.set_title("Time Mean Zonal Velocity (m/s)")
uVals = np.arange(-0.1,0.3+0.01,0.02)
uValsCbar = np.arange(-0.1,0.31,0.1)
print densityLayers

thicknessAveU = np.zeros((uTimeMean[1,...]).shape,dtype='float64')
for iLayer in range(0,zl.shape[0]):
	thicknessAveU[:,:] = thicknessAveU[:,:] + (2.0*refThickness[2.0*iLayer] * uTimeMean[iLayer,:,:])
thicknessAveU = thicknessAveU / refThickness[::2].sum()

cs = ax.contourf(x,y,thicknessAveU[:,:],uVals,cmap=plt.cm.jet)
ax.annotate('a',xy=(0.01,1.025),xycoords='axes fraction',fontsize=15)


cbar = fig.colorbar(cs,ax=ax, ticks=uValsCbar)
cbar.ax.set_title('m/s')

ax.contour(x,y,topo[1::,:],10,colors='k')
#ax.plot(x,y[upperIndex]*np.ones(x.shape[0],dtype='float64'),'k')
#ax.plot(x,y[lowerIndex]*np.ones(x.shape[0],dtype='float64'),'k')
ax.set_ylabel('y (km)')
ax.get_xaxis().set_visible(False)


ax  = fig.add_subplot(3,1,2)
vVals = np.arange(-0.1,0.1+0.01,0.01)
vValsCbar = np.arange(-0.1,0.11,0.05)

thicknessAveV = np.zeros((vTimeMean[1,...]).shape,dtype='float64')
for iLayer in range(0,zl.shape[0]):
	thicknessAveV[:,:] = thicknessAveV[:,:] + 2.0*(refThickness[2.0*iLayer] * vTimeMean[iLayer,:,:])
thicknessAveV = thicknessAveV / refThickness[::2].sum()
cs = ax.contourf(x,y,thicknessAveV,15,cmap=plt.cm.RdBu_r)
cbar = fig.colorbar(cs,ax=ax)#,ticks = vValsCbar)
cbar.ax.set_title('m/s')
ax.annotate('b',xy=(0.01,1.025),xycoords='axes fraction',fontsize=15)

ax.contour(x,y,topo[1::,:],10,colors='k')

#ax.plot(x,y[upperIndex]*np.ones(x.shape[0],dtype='float64'),'k')
#ax.plot(x,y[lowerIndex]*np.ones(x.shape[0],dtype='float64'),'k')
ax.set_ylabel('y (km)')
ax.get_xaxis().set_visible(False)



ax  = fig.add_subplot(3,1,3)
ekeVals = np.arange(0,50.5,2) 
ekeValsCbar = np.arange(0,51,10)
thicknessAveEKE = np.zeros((ekeTimeMean[1,...]).shape,dtype='float64')
for iLayer in range(0,zl.shape[0]):
	thicknessAveEKE[:,:] = thicknessAveEKE[:,:] + (2.0*densityLayers[2.0*iLayer]*refThickness[2.0*iLayer] * ekeTimeMean[iLayer,:,:])
thicknessAveEKE = thicknessAveEKE / refThickness[::2].sum()

cs = ax.contourf(x,y,densityLayers[6]*ekeTimeMean[2,:,:],15,cmap=plt.cm.jet)
ax.annotate('c',xy=(0.01,1.025),xycoords='axes fraction',fontsize=15)


cbar = fig.colorbar(cs,ax=ax)#,ticks=ekeValsCbar)
cbar.ax.set_title('J.m$^{-3}$')

ax.contour(x,y,topo[1::,:],10,colors='k')

#ax.plot(x,y[upperIndex]*np.ones(x.shape[0],dtype='float64'),'k')
#ax.plot(x,y[lowerIndex]*np.ones(x.shape[0],dtype='float64'),'k')
ax.set_xlabel('x (km)')
ax.set_ylabel('y (km)')


#plt.show()


transientWaveActivityFileName = 'waveActivity_transient_RandTopoSmooth.nc'
waveFile = scipy_netcdf.netcdf_file(outputDirectory + transientWaveActivityFileName, 'r')
EPVector = []
EPVector.append(waveFile.variables['EP_X'][:,:,:])
EPVector.append(waveFile.variables['EP_Y'][:,:,:])
EPVector.append(waveFile.variables['EP_Z'][:,:,:])

waveFile.close()

normEP = np.sqrt(EPVector[0]**2+EPVector[1]**2)
print normEP.mean
threshold = 0.0001
maskArray = normEP[iLayer,:,:] < threshold
EP_X_mask = np.ma.masked_array(EPVector[0][3,:,:], mask=maskArray) 
EP_Y_mask = np.ma.masked_array(EPVector[1][3,:,:], mask=maskArray) 
#div4 = computeDivergence(EPVector[0],EPVector[1], deltaX*1.0e3, deltaY*1.0e3, layerToPlot)

fig = plt.figure(2)
plotInterval = 20 


ekeVals = np.arange(0,50.5,2) 
ekeValsCbar = np.arange(0,51,10)
ax  = fig.add_subplot(1,1,1, aspect='equal')
cs = ax.contourf(x,y[1::],densityLayers[4]*ekeTimeMean[2,1:ekeTimeMean.shape[1],:],15,cmap=plt.cm.jet)
cbar_ax = fig.add_axes([0.91, 0.24, 0.015, 0.52])
cbar = fig.colorbar(cs, cax=cbar_ax)##,ticks=cbarTicks)
cbar.ax.set_title(r'(J.m$^{-3}$)')


#cbar = fig.colorbar(cs,ax=ax)#,ticks=ekeValsCbar)
#$cbar.ax.set_title('J.m$^{-3}$')
#ax.contour(x,y,topo6[1:topo6.shape[0]-1,:],15,colors='k')
X,Y = np.meshgrid(x,y)


Q = ax.quiver( X[::plotInterval, ::plotInterval], Y[::plotInterval, ::plotInterval], 1.0*EP_X_mask[::plotInterval, ::plotInterval], 
                                          EP_Y_mask[::plotInterval, ::plotInterval], 
										  pivot='mid', color='k',units='x',scale=0.3e-4,scale_units='x')
qk = plt.quiverkey(Q, 0.02, 1.01, 3.0*threshold*25.0, r'10J.m$^{-3}$', coordinates='axes',fontproperties={'weight': 'bold'})
ax.set_xlabel('x (km)')
ax.set_ylabel('y (km)')
print "Threshold", 3.0*threshold





transientFluxFileName = 'transient_fluxes_randTopoSmooth.nc'

transientFluxFile = scipy_netcdf.netcdf_file(outputDirectory + transientFluxFileName, 'r')

uFluxTransient  = transientFluxFile.variables['uEddyFlux'][:,:,:]
vFluxTransient  = transientFluxFile.variables['vEddyFlux'][:,:,:]

uFluxTotal  = transientFluxFile.variables['uTotalFlux'][:,:,:]
vFluxTotal  = transientFluxFile.variables['vTotalFlux'][:,:,:]


uFluxStanding = uFluxTotal - uFluxTransient
vFluxStanding = vFluxTotal - vFluxTransient
vTimeMean = transientFluxFile.variables['v_time_mean'][:,:,:]



uFluxGradStanding = np.gradient(uFluxStanding,1.0,deltaY*1.0e3,deltaX*1.0e3)
vFluxGradStanding = np.gradient(vFluxStanding,1.0,deltaY*1.0e3,deltaX*1.0e3)


uFluxGradTransient = np.gradient(uFluxTransient,1.0,deltaY*1.0e3,deltaX*1.0e3)
vFluxGradTransient = np.gradient(vFluxTransient,1.0,deltaY*1.0e3,deltaX*1.0e3)

uFluxGradTotal = np.gradient(uFluxTotal,1.0,deltaY*1.0e3,deltaX*1.0e3)
vFluxGradTotal = np.gradient(vFluxTotal,1.0,deltaY*1.0e3,deltaX*1.0e3)




divStandingFlux  = uFluxGradStanding[2] + vFluxGradStanding[1]
divTransientFlux = uFluxGradTransient[2] + vFluxGradTransient[1]
divTotalFlux = uFluxGradTotal[2] + vFluxGradTotal[1]

print vFluxTransient[iLayer,:,:].shape
print y.shape

fig = plt.figure(3)
ax  = fig.add_subplot(1,1,1)
ax.plot(y,vFluxTotal[iLayer,:,:].mean(axis=-1),'k')
ax.plot(y,vFluxTransient[iLayer,:,:].mean(axis=-1),'b')
ax.plot(y,(vFluxStanding[iLayer,:,:]).mean(axis=-1),'r')


'''
r = 0.75
startNumber = -0.12
plotLevels = [startNumber]
nPoints =10
for i in range(0,nPoints):
	plotLevels.append(plotLevels[-1] * r)
plotLevels = np.asarray(plotLevels + [0.0] + plotLevels[::-1])
plotLevels[nPoints+2::] = -1.0*plotLevels[nPoints+2::] 
print plotLevels

r = 0.5
cbarTicks = [startNumber]
nPoints =3
for i in range(0,nPoints):
	cbarTicks.append(cbarTicks[-1] * r)
cbarTicks = np.asarray(cbarTicks + [0.0] + cbarTicks[::-1])
cbarTicks[nPoints+2::] = -1.0*cbarTicks[nPoints+2::] 
'''






fig = plt.figure(4)
ax  = fig.add_subplot(3,1,1)
cs = ax.contourf(x,y,(divStandingFlux[iLayer,:,:])*1.0e3,15,cmap=plt.cm.jet)
#ax.contour(x,y,topo[1::,...],15,colors='k')
fig.colorbar(cs, ax=ax)

ax.annotate('a',xy=(0.01,1.025),xycoords='axes fraction',fontsize=15)
#ax.get_xaxis().set_visible(False)

print y.shape
print divStandingFlux.shape
print divTransientFlux.shape

ax  = fig.add_subplot(3,1,2)
cs = ax.contourf(x,y,(divTransientFlux[iLayer,:,:])*1.0e3,15,cmap=plt.cm.jet)
fig.colorbar(cs, ax=ax)

#ax.contour(x,y,topo[1::,:],15,colors='k')
ax.annotate('b',xy=(0.01,1.025),xycoords='axes fraction',fontsize=15)
ax.set_ylabel('y (km)')
#ax.get_xaxis().set_visible(False)



ax  = fig.add_subplot(3,1,3)
cs = ax.contourf(x,y,(divTotalFlux[iLayer,:,:])*1.0e3,15,cmap=plt.cm.jet)
fig.colorbar(cs, ax=ax)

#ax.contour(x,y,topo[1::,:],15,colors='k')
ax.annotate('c',xy=(0.01,1.025),xycoords='axes fraction',fontsize=15)
ax.set_xlabel('x (km)')


#fig.subplots_adjust(right=0.8)
#cbar_ax = fig.add_axes([0.85, 0.08, 0.05, 0.82])
#fig.colorbar(cs, cax=cbar_ax,ticks=cbarTicks)
#cbar_ax.set_title(r'$\nabla \cdot (\vec{u}z)$ (m.s$^{-1}$) ')


zonalWaveActivityFileName = 'waveActivity_Zonal_randTopoSmooth1.nc'

waveFile = scipy_netcdf.netcdf_file(outputDirectory+zonalWaveActivityFileName, 'r')

EPVector_Z = []
EPVector_Z.append(waveFile.variables['EP_X'][:,:,:])
EPVector_Z.append(waveFile.variables['EP_Y'][:,:,:])
EPVector_Z.append(waveFile.variables['EP_Z'][:,:,:])
print EPVector_Z[2][iLayer,:,:].shape
print x[1:-1].shape
print y.shape


iLayer = 2
fig = plt.figure(6)
#ax  = fig.add_subplot(3,1,1)
#cs = ax.contourf(x,y,ekeTimeMean[iLayer,:,:]* densityLayers[iLayer],15,cmap=plt.cm.RdBu_r)
#fig.colorbar(cs,ax=ax)

lowerLatitude = 500.0
upperLatitude = 1200.0
lowerIndex = np.nonzero(y >= lowerLatitude)[0][0]
upperIndex = np.nonzero(y >= upperLatitude)[0][0]

thicknessAveW = np.zeros(EPVector_Z[2][1,:,:].shape,dtype='float64')

for iLayer in range(0,EPVector_Z[2].shape[0]):
	thicknessAveW[:,:] = thicknessAveW[:,:] + (densityLayers[iLayer]*refThickness[iLayer] * EPVector_Z[2][iLayer,:,:])
thicknessAveW = thicknessAveW / refThickness.sum()

ax  = fig.add_subplot(2,1,1)
ln1 = ax.plot(x[1:-1],-EPVector_Z[2][3,lowerIndex:upperIndex,:].mean(axis=0),'k-',label = 'Wave Flux')
ax.set_ylabel(r'$W_{\rho}$ (J.m$^{-3}$)')
#yfm = ax.yaxis.get_major_formatter()
#yfm.set_powerlimits([ -2, 0])
ax.annotate('a',xy=(0.98,1.05),xycoords='axes fraction',fontsize=15)

ax2 = ax.twinx()
ln2 = ax2.plot(x[1:-1],vTimeMean[2,lowerIndex:upperIndex,1:-1].mean(axis=0),'k-.',label='V')
ax2.set_ylabel('V (m/s)')
ax.get_xaxis().set_visible(False)
ax.set_xlim(x[0],x[-1])
# added these three lines
lns = ln1+ln2
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc='lower left', shadow=True, prop={'size':10})


ax3  = fig.add_subplot(2,1,2)
ln3 = ax3.plot(x[1:-1], densityLayers[4]*ekeTimeMean[2,lowerIndex:upperIndex,1:-1].mean(axis=0),'k',label='EKE')
ax3.set_xlabel('x (km)')
ax3.set_xlim(x[0],x[-1])
ax3.annotate('b',xy=(0.98,1.025),xycoords='axes fraction',fontsize=15)
ax3.set_ylabel('EKE (J.m$^{-3}$)')
ax4 = ax3.twinx()
ln4 = ax4.plot(x[1:-1], topo[lowerIndex:upperIndex,1:-1].mean(axis=0),'k--',label='Topo')
lns = ln3+ln4

labs = [l.get_label() for l in lns]
ax3.legend(lns, labs, loc='upper left', shadow=True, prop={'size':10})
ax4.set_yticks([1000,1500,2000])
ax4.set_ylabel('Topo height (m)')
ax3.set_xlim(x[0],x[-1])
ax4.set_xlim(x[0],x[-1])




fig = plt.figure(7)
ax  = fig.add_subplot(1,1,1)
cs = ax.contourf(x[1:-1],y[1::],EPVector_Z[2][iLayer,:,:],15,cmap=plt.cm.jet)
fig.colorbar(cs,ax=ax)

cs = ax.contour(x,y,vTimeMean[iLayer,:,:],15,colors='k')
ax.set_xlim(x[0],x[-1])



#vValues = np.arange(-0.2,0.15+0.01,0.05)
#ax.contour(x,y,vTimeMean[iLayer,:,:],vValues,colors='k')
#ax.contour(x,y,topo5,15,colors='k',linestyles='dotted')
#cbar = fig.colorbar(cs,ax=ax,ticks=cbarTicks)
ax.set_ylabel('y (km)')
ax.set_xlabel('x (km)')
#cbar.ax.set_title('EP Flux (kJ.m$^{-2}$)')





waveFile.close()

fig = plt.figure(8)

#ax  = fig.add_subplot(3,1,1)
#cs = ax.contourf(x,y,ekeTimeMean[iLayer,:,:]* densityLayers[iLayer],15,cmap=plt.cm.RdBu_r)
#fig.colorbar(cs,ax=ax)

topoVals = np.arange(0.0,3005.0,200.0)
topoTicks = np.arange(0.0,3100.0,1000.0)

ax  = fig.add_subplot(1,1,1, aspect='equal')
cs = ax.contourf(x,y,(topo[1::,...]),topoVals,cmap=plt.cm.jet)
cbar = fig.colorbar(cs,ax=ax,ticks=topoTicks)
cbar.ax.set_title('m')

ax.set_ylabel('y (km)')
ax.set_xlabel('x (km)')




plt.show()

