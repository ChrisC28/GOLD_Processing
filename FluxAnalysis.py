import numpy as np
import os
from scipy.io import netcdf as scipy_netcdf
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm

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
outputDirectory = '/home/chris/GOLD/IsoTopo_4_Processing/'

standingFluxFileName  = 'standing_fluxes_isoTopo4.nc'
transientFluxFileName = 'transient_fluxes_isoTopo4.nc'
EKE_FileName = 'postprocessed_IsoTopo4.nc'

topoFileName   = 'isoTopo_4.nc'

standingFluxFile  = scipy_netcdf.netcdf_file(outputDirectory + standingFluxFileName, 'r')
transientFluxFile = scipy_netcdf.netcdf_file(outputDirectory + transientFluxFileName, 'r')
ekeFile           = scipy_netcdf.netcdf_file(outputDirectory + EKE_FileName, 'r')

topoFile  = scipy_netcdf.netcdf_file(outputDirectory + topoFileName, 'r')
#Get dimensions
xName = 'x'
yName = 'y'
zlName = 'zl'
ziName = 'zi'

x  = standingFluxFile.variables[xName][:]
y  = standingFluxFile.variables[yName][:]
zl = standingFluxFile.variables[zlName][:]
nZ = zl.shape[0]

deltaX = x[1]-x[0]
deltaY = y[1]-y[0]


#uFluxStanding   = standingFluxFile.variables['statEddyFluxU'][:,:,:]
#vFluxStanding   = standingFluxFile.variables['statEddyFluxV'][:,:,:]

uFluxTransient  = transientFluxFile.variables['uEddyFlux'][:,:,:]
vFluxTransient  = transientFluxFile.variables['vEddyFlux'][:,:,:]

uFluxTotal  = transientFluxFile.variables['uTotalFlux'][:,:,:]
vFluxTotal  = transientFluxFile.variables['vTotalFlux'][:,:,:]


uFluxStanding = uFluxTotal - uFluxTransient
vFluxStanding = vFluxTotal - vFluxTransient

topo        = meanDepth - topoFile.variables['topo'][:,:]

ekeTimeMean  = ekeFile.variables['eke'][:,:,:]



ekeFile.close()
standingFluxFile.close()
transientFluxFile.close()
topoFile.close()

uFluxGradStanding = np.gradient(uFluxStanding,1.0,deltaY*1.0e3,deltaX*1.0e3)
vFluxGradStanding = np.gradient(vFluxStanding,1.0,deltaY*1.0e3,deltaX*1.0e3)


uFluxGradTransient = np.gradient(uFluxTransient,1.0,deltaY*1.0e3,deltaX*1.0e3)
vFluxGradTransient = np.gradient(vFluxTransient,1.0,deltaY*1.0e3,deltaX*1.0e3)

uFluxGradTotal = np.gradient(uFluxTotal,1.0,deltaY*1.0e3,deltaX*1.0e3)
vFluxGradTotal = np.gradient(vFluxTotal,1.0,deltaY*1.0e3,deltaX*1.0e3)


lowerLatitude  = 700.0
upperLatitude = 1300.0
lowerIndex = np.nonzero(y >= lowerLatitude)[0][0]
upperIndex = np.nonzero(y >= upperLatitude)[0][0]

yAveragedEKE = np.mean(ekeTimeMean[:,lowerIndex:upperIndex,:],axis=1)
averagedEKE = np.zeros(x.shape[0],dtype='float64')
for iZ in range(0,nZ/2):
	averagedEKE[:] = averagedEKE[:] + 2.0*(yAveragedEKE[iZ,:] * densityLayers[2.0*iZ])






divStandingFlux  = uFluxGradStanding[2] + vFluxGradStanding[1]
divTransientFlux = uFluxGradTransient[2] + vFluxGradTransient[1]
divTotalFlux = uFluxGradTotal[2] + vFluxGradTotal[1]






averagedStandingDiv   = np.mean(divStandingFlux[:,lowerIndex:upperIndex,:],axis=1)
averagedTransientDiv  = np.mean(divTransientFlux[:,lowerIndex:upperIndex,:],axis=1)


totalAveStandingDiv= np.zeros(x.shape[0],dtype='float64')
totalAveTransientDiv= np.zeros(x.shape[0],dtype='float64')


for iZ in range(0,nZ):
	totalAveStandingDiv[:] = totalAveStandingDiv[:] + (averagedStandingDiv[iZ,...] * densityLayers[iZ])
	totalAveTransientDiv[:] = totalAveTransientDiv[:] + (averagedTransientDiv[iZ,...] * densityLayers[iZ])
	
print averagedEKE.shape
print x.shape


vertIntDivStanding  =   np.zeros([divStandingFlux.shape[1],divStandingFlux.shape[2]]  ,dtype='float64')
vertIntDivTransient =   np.zeros([divTransientFlux.shape[1],divTransientFlux.shape[2]],dtype='float64')
print vertIntDivTransient.shape
for iLayer in range(0,averagedStandingDiv.shape[0]):
	vertIntDivStanding[:,:]  = vertIntDivStanding[:,:] + (divStandingFlux[iLayer,...] * refThickness[iLayer])
	vertIntDivTransient[:,:] = vertIntDivTransient[:,:] + (divTransientFlux[iLayer,...] * refThickness[iLayer])










divStandingFluxYAve = np.mean(np.gradient(divStandingFlux[:,lowerIndex:upperIndex,:],1.0,deltaY*1.0e3,deltaX*1.0e3)[1],axis=1)
divTransientFluxYAve = np.mean(np.gradient(divTransientFlux[:,lowerIndex:upperIndex,:],1.0,deltaY*1.0e3,deltaX*1.0e3)[1],axis=1)


layerToPlot = 3
fig = plt.figure(0)
ax1  = fig.add_subplot(1,1,1)
ax1.plot(x,yAveragedEKE[layerToPlot/2,...]*deltaY*1.0e3,'k')

ax2 = ax1.twinx()
yAve = divStandingFlux[:]
ax2.plot(x,divStandingFluxYAve[layerToPlot,...]*deltaY*1.0e3,'b')
ax2.plot(x,divTransientFluxYAve[layerToPlot,...]*deltaY*1.0e3,'r')








#plotLevels = [-0.12,-0.09,-0.06,-0.045,-0.03,-0.015,-7.5e-3,0.0,7.5e-3,0.015,0.03,0.06,0.12]
#plotLevels = np.logspace(-0.12, 0.12, num=50, endpoint=True, base=10.0
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
	



fig = plt.figure(1)
ax  = fig.add_subplot(3,1,1)
cs = ax.contourf(x,y,(divStandingFlux[layerToPlot,1::,:])*1.0e3,plotLevels,cmap=plt.cm.jet)
ax.contour(x,y,topo[1:-1,...],15,colors='k')
ax.annotate('a',xy=(0.01,1.025),xycoords='axes fraction',fontsize=15)
ax.get_xaxis().set_visible(False)

print y.shape
print divStandingFlux.shape
print divTransientFlux.shape

ax  = fig.add_subplot(3,1,2)
cs = ax.contourf(x,y,(divTransientFlux[layerToPlot,1::,:])*1.0e3,plotLevels,cmap=plt.cm.jet)
ax.contour(x,y,topo[1:-1,:],15,colors='k')
ax.annotate('b',xy=(0.01,1.025),xycoords='axes fraction',fontsize=15)
ax.set_ylabel('y (km)')
ax.get_xaxis().set_visible(False)



ax  = fig.add_subplot(3,1,3)
cs = ax.contourf(x,y,(divTotalFlux[layerToPlot,1::,:])*1.0e3,plotLevels,cmap=plt.cm.jet)
ax.contour(x,y,topo[1:-1,:],15,colors='k')
ax.annotate('c',xy=(0.01,1.025),xycoords='axes fraction',fontsize=15)
ax.set_xlabel('x (km)')

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.08, 0.05, 0.82])
fig.colorbar(cs, cax=cbar_ax,ticks=cbarTicks)
cbar_ax.set_title(r'$\nabla \cdot (\vec{u}z)$ (m.s$^{-1}$) ')






fig = plt.figure(4)
ax  = fig.add_subplot(1,1,1)
ax.plot(y,vFluxTotal[layerToPlot,1::,:].mean(axis=-1),'k')
ax.plot(y,vFluxTransient[layerToPlot,1::,:].mean(axis=-1),'b')
ax.plot(y,(vFluxStanding[layerToPlot,1::,:]).mean(axis=-1),'r')


fig = plt.figure(5)
ax  = fig.add_subplot(1,1,1)
cs = ax.contourf(x,y,vFluxStanding[layerToPlot,1::,:],15,cmap=plt.cm.jet)
fig.colorbar(cs,ax=ax)
ax.contour(x,y,topo[1:-1,...],15,colors='k')

print y.shape
print divStandingFlux.shape
print divTransientFlux.shape
fig = plt.figure(6)
ax  = fig.add_subplot(1,1,1)
cs = ax.contourf(x,y,vFluxTransient[layerToPlot,1::,:],15,cmap=plt.cm.jet)
fig.colorbar(cs,ax=ax)
ax.contour(x,y,topo[1:-1,:],15,colors='k')

fig = plt.figure(7)
ax  = fig.add_subplot(1,1,1)
cs = ax.contourf(x,y,vFluxTotal[layerToPlot,1::,:],15,cmap=plt.cm.jet)
fig.colorbar(cs,ax=ax)
ax.contour(x,y,topo[1:-1,:],15,colors='k')


plt.show()
