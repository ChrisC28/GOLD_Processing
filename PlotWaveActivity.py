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
dirName = '/home/chris/GOLD/IsoTopo_' + str(iDir) + '_Processing/'
zonalWaveActivityFileName = dirName + 'waveActivity_Transient_isoTopo' + str(iDir) + '.nc'
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

threshold = 0.01
normEP5 = np.sqrt(EPVector_5[0]**2+EPVector_5[1]**2)
maskArray = normEP5[layerToPlot,:,:] < threshold
EP_X_mask5 = np.ma.masked_array(EPVector_5[0][layerToPlot,:,:], mask=maskArray) 
EP_Y_mask5 = np.ma.masked_array(EPVector_5[1][layerToPlot,:,:], mask=maskArray) 

div5 = computeDivergence(EPVector_5[0],EPVector_5[1], deltaX*1.0e3, deltaY*1.0e3, layerToPlot)

###########################################################
iDir = 6
dirName = '/home/chris/GOLD/IsoTopo_' + str(iDir) + '_Processing/'
zonalWaveActivityFileName = dirName + 'waveActivityTransient_isoTopo' + str(iDir) + '.nc'
topoFileName  = dirName + 'isoTopo_' + str(iDir) + '.nc'
EKE_fileName = dirName + 'postprocessed_IsoTopo' + str(iDir) + '.nc'
print EKE_fileName
	

waveFile = scipy_netcdf.netcdf_file(zonalWaveActivityFileName, 'r')
topoFile  = scipy_netcdf.netcdf_file(topoFileName, 'r')
ekeFile  = scipy_netcdf.netcdf_file(EKE_fileName, 'r')

topo6        = meanDepth - topoFile.variables[topoName][:,:]
topoFile.close()


EPVector_6 = []


EPVector_6.append(waveFile.variables['EP_X'][:,:,:])
EPVector_6.append(waveFile.variables['EP_Y'][:,:,:])
EPVector_6.append(waveFile.variables['EP_Z'][:,:,:])

ekeTimeMean6 = ekeFile.variables['eke'][:,:,:]
vTimeMean6 = ekeFile.variables['v_time_mean'][:,:,:]

ekeFile.close()
waveFile.close()

normEP6 = np.sqrt(EPVector_6[0]**2+EPVector_6[1]**2)

threshold6 = 0.075
maskArray = normEP6[layerToPlot,:,:] < threshold
EP_X_mask6 = np.ma.masked_array(EPVector_6[0][layerToPlot,:,:], mask=maskArray) 
EP_Y_mask6 = np.ma.masked_array(EPVector_6[1][layerToPlot,:,:], mask=maskArray) 
div6 = computeDivergence(EPVector_6[0],EPVector_6[1], deltaX*1.0e3, deltaY*1.0e3, layerToPlot)


###########################################################
iDir = 7
dirName = '/home/chris/GOLD/IsoTopo_' + str(iDir) + '_Processing/'
zonalWaveActivityFileName = dirName + 'waveActivity_transient_isoTopo' + str(iDir) + '.nc'
topoFileName  = dirName + 'isoTopo_' + str(iDir) + '.nc'
EKE_fileName = dirName + 'postprocessed_IsoTopo' + str(iDir) + '.nc'
print EKE_fileName
	

waveFile = scipy_netcdf.netcdf_file(zonalWaveActivityFileName, 'r')
topoFile  = scipy_netcdf.netcdf_file(topoFileName, 'r')
ekeFile  = scipy_netcdf.netcdf_file(EKE_fileName, 'r')

topo7        = meanDepth - topoFile.variables[topoName][:,:]
topoFile.close()


EPVector_7 = []


EPVector_7.append(waveFile.variables['EP_X'][:,:,:])
EPVector_7.append(waveFile.variables['EP_Y'][:,:,:])
EPVector_7.append(waveFile.variables['EP_Z'][:,:,:])

ekeTimeMean7 = ekeFile.variables['eke'][:,:,:]
vTimeMean7 = ekeFile.variables['v_time_mean'][:,:,:]

ekeFile.close()
waveFile.close()

normEP7 = np.sqrt(EPVector_7[0]**2+EPVector_7[1]**2)

maskArray = normEP7[layerToPlot,:,:] < threshold
EP_X_mask7 = np.ma.masked_array(EPVector_7[0][layerToPlot,:,:], mask=maskArray)
EP_Y_mask7 = np.ma.masked_array(EPVector_7[1][layerToPlot,:,:], mask=maskArray) 
div7 = computeDivergence(EPVector_7[0],EPVector_7[1], deltaX*1.0e3, deltaY*1.0e3, layerToPlot)

#####################################################################

iDir = 8
dirName = '/home/chris/GOLD/IsoTopo_' + str(iDir) + '_Processing/'
zonalWaveActivityFileName = dirName + 'waveActivity_Transient_isoTopo' + str(iDir) + '.nc'
topoFileName  = dirName + 'isoTopo_' + str(iDir) + '.nc'
EKE_fileName = dirName + 'postprocessed_IsoTopo' + str(iDir) + '.nc'
print EKE_fileName
	

waveFile = scipy_netcdf.netcdf_file(zonalWaveActivityFileName, 'r')
topoFile  = scipy_netcdf.netcdf_file(topoFileName, 'r')
ekeFile  = scipy_netcdf.netcdf_file(EKE_fileName, 'r')

topo8        = meanDepth - topoFile.variables[topoName][:,:]
topoFile.close()


EPVector_8 = []


EPVector_8.append(waveFile.variables['EP_X'][:,:,:])
EPVector_8.append(waveFile.variables['EP_Y'][:,:,:])
EPVector_8.append(waveFile.variables['EP_Z'][:,:,:])

ekeTimeMean8 = ekeFile.variables['eke'][:,:,:]
vTimeMean8 = ekeFile.variables['v_time_mean'][:,:,:]

ekeFile.close()
waveFile.close()

normEP8 = np.sqrt(EPVector_8[0]**2+EPVector_7[1]**2)

maskArray = normEP8[layerToPlot,:,:] < threshold
EP_X_mask8 = np.ma.masked_array(EPVector_7[0][layerToPlot,:,:], mask=maskArray)
EP_Y_mask8 = np.ma.masked_array(EPVector_7[1][layerToPlot,:,:], mask=maskArray)


div8 = computeDivergence(EPVector_8[0],EPVector_8[1], deltaX*1.0e3, deltaY*1.0e3, layerToPlot)




#########################################################

iDir = 4
dirName = '/home/chris/GOLD/IsoTopo_' + str(iDir) + '_Processing/'
zonalWaveActivityFileName = dirName + 'waveActivityTransient_isoTopo' + str(iDir) + '.nc'
topoFileName  = dirName + 'isoTopo_' + str(iDir) + '.nc'
EKE_fileName = dirName + 'postprocessed_IsoTopo' + str(iDir) + '.nc'
print EKE_fileName
	
print zonalWaveActivityFileName
waveFile = scipy_netcdf.netcdf_file(zonalWaveActivityFileName, 'r')
topoFile  = scipy_netcdf.netcdf_file(topoFileName, 'r')
ekeFile  = scipy_netcdf.netcdf_file(EKE_fileName, 'r')

topo4        = meanDepth - topoFile.variables[topoName][:,:]
topoFile.close()


EPVector_4 = []


EPVector_4.append(waveFile.variables['EP_X'][:,:,:])
EPVector_4.append(waveFile.variables['EP_Y'][:,:,:])
EPVector_4.append(waveFile.variables['EP_Z'][:,:,:])

ekeTimeMean4 = ekeFile.variables['eke'][:,:,:]
vTimeMean4 = ekeFile.variables['v_time_mean'][:,:,:]

ekeFile.close()
waveFile.close()

normEP4 = np.sqrt(EPVector_6[0]**2+EPVector_6[1]**2)

threshold4 = 0.01
maskArray = normEP4[layerToPlot,:,:] < threshold
EP_X_mask4 = np.ma.masked_array(EPVector_4[0][layerToPlot,:,:], mask=maskArray) 
EP_Y_mask4 = np.ma.masked_array(EPVector_4[1][layerToPlot,:,:], mask=maskArray) 
div4 = computeDivergence(EPVector_4[0],EPVector_4[1], deltaX*1.0e3, deltaY*1.0e3, layerToPlot)

#########################################################




print x.shape
print y.shape
print vTimeMean6[layerToPlot,:,:].shape
print topo6[1::,:].shape
print ekeTimeMean6[layerToPlot,1:ekeTimeMean6.shape[1],:].shape
X,Y = np.meshgrid(x,y)

fig = plt.figure(1)
ax  = fig.add_subplot(3,1,1)
cs  = ax.contourf(x,y,EPVector_6[2][layerToPlot,...],15,cmap=plt.cm.jet)
fig.colorbar(cs,ax=ax)

#ax.contour(x,y,vTimeMean6[layerToPlot,1::,1:vTimeMean6.shape[2]-1],15,colors='k')
ax.contour(x,y,vTimeMean6[layerToPlot,1::,0:vTimeMean6.shape[2]],15,colors='k')
#ax.contour(x,y,topo6[1:topo6.shape[0]-1,1:topo6.shape[1]-1],15,colors='k')
ax.contour(x,y,topo6[1:topo6.shape[0]-1,:],15,colors='k')




ax  = fig.add_subplot(3,1,2)
ax.contourf(x,y,EPVector_5[2][layerToPlot,...],15,cmap=plt.cm.jet)
ax.contour(x,y,vTimeMean5[layerToPlot,1::,0:vTimeMean6.shape[2]],15,colors='k')
ax.contour(x,y,topo5[1:topo6.shape[0]-1,:],15,colors='k')
fig.colorbar(cs,ax=ax)


ax  = fig.add_subplot(3,1,3)
ax.contourf(x,y,EPVector_7[2][layerToPlot,...],15,cmap=plt.cm.jet)
ax.contour(x,y,vTimeMean7[layerToPlot,1::,0:vTimeMean6.shape[2]],15,colors='k')
ax.contour(x,y,topo7[1:topo6.shape[0]-1,:],15,colors='k')
fig.colorbar(cs,ax=ax)


fig = plt.figure(2)
ax  = fig.add_subplot(3,1,1)
ax.contourf(x,y,EPVector_6[2][layerToPlot,...],15,cmap=plt.cm.jet)
ax.contour(x,y,ekeTimeMean6[layerToPlot,1:ekeTimeMean6.shape[1],:],15,colors='k')

ax  = fig.add_subplot(3,1,2)
ax.contourf(x,y,EPVector_5[2][layerToPlot,...],15,cmap=plt.cm.jet)
ax.contour(x,y,ekeTimeMean5[layerToPlot,1:ekeTimeMean6.shape[1],:],15,colors='k')

ax  = fig.add_subplot(3,1,3)
ax.contourf(x,y,EPVector_7[2][layerToPlot,...],15,cmap=plt.cm.jet)
ax.contour(x,y,ekeTimeMean7[layerToPlot,1:ekeTimeMean6.shape[1],:],15,colors='k')

fig = plt.figure(3)
plotInterval = 20 


cbarTicksToPlot = np.arange(0,18.1,1.0)
cbarTicks       = np.arange(0,18.1,3.0)
thicknessAveEKE6 = np.zeros((ekeTimeMean6[1,...]).shape,dtype='float64')
for iLayer in range(0,zl.shape[0]/2):
	thicknessAveEKE6[:,:] = thicknessAveEKE6[:,:] + (densityLayers[2.0*iLayer]*refThickness[2.0*iLayer] * ekeTimeMean6[iLayer,:,:])
thicknessAveEKE6 = thicknessAveEKE6 / refThickness[::2].sum()

ax  = fig.add_subplot(3,1,1)
cs = ax.contourf(x,y,thicknessAveEKE6[1:ekeTimeMean6.shape[1],:],15,cmap=plt.cm.RdBu_r)
#fig.colorbar(cs,ax=ax,ticks=cbarTicks)
ax.contour(x,y,topo6[1:topo6.shape[0]-1,:],15,colors='k')



Q = ax.quiver( X[::plotInterval, ::plotInterval], Y[::plotInterval, ::plotInterval], 1.0*EP_X_mask6[::plotInterval, ::plotInterval], 
                                          2.0*EP_Y_mask6[::plotInterval, ::plotInterval], 
										  pivot='mid', color='k',units='x',scale=0.15e-3,scale_units='x')
qk = plt.quiverkey(Q, 0.02, 1.1, 3.0*threshold, r'30J.m$^{-3}$', coordinates='axes',fontproperties={'weight': 'bold'})
print normEP6[layerToPlot,:,:].std()
print normEP5[layerToPlot,:,:].std()
print normEP7[layerToPlot,:,:].std()
ax.get_xaxis().set_visible(False)


ax  = fig.add_subplot(3,1,2)
thicknessAveEKE5 = np.zeros((ekeTimeMean5[1,...]).shape,dtype='float64')
for iLayer in range(0,zl.shape[0]/2):
	thicknessAveEKE5[:,:] = thicknessAveEKE5[:,:] + (densityLayers[2.0*iLayer]*refThickness[2.0*iLayer] * ekeTimeMean5[iLayer,:,:])
thicknessAveEKE5 = thicknessAveEKE5 / refThickness[::2].sum()

ax.contourf(x,y,thicknessAveEKE5[1:ekeTimeMean6.shape[1],:],15,cmap=plt.cm.RdBu_r)
#fig.colorbar(cs,ax=ax,ticks=cbarTicks)

ax.contour(x,y,topo5[1:topo6.shape[0]-1,:],15,colors='k')
Q = ax.quiver( X[::plotInterval, ::plotInterval], Y[::plotInterval, ::plotInterval], EP_X_mask5[::plotInterval, ::plotInterval], 
                                          EP_Y_mask5[::plotInterval, ::plotInterval], 
										  pivot='mid', color='k',units='x',scale=0.15e-3,scale_units='x')



#qk = plt.quiverkey(Q, 0.02, 1.1, 3.0*threshold, r'30', coordinates='axes',fontproperties={'weight': 'bold'})
ax.get_xaxis().set_visible(False)
ax.set_ylabel('y (km)')


ax  = fig.add_subplot(3,1,3)
thicknessAveEKE7 = np.zeros((ekeTimeMean7[1,...]).shape,dtype='float64')
for iLayer in range(0,zl.shape[0]/2):
	thicknessAveEKE7[:,:] = thicknessAveEKE7[:,:] + (densityLayers[2.0*iLayer]*refThickness[2.0*iLayer] * ekeTimeMean7[iLayer,:,:])
thicknessAveEKE7 = thicknessAveEKE7 / refThickness[::2].sum()

cs = ax.contourf(x,y,thicknessAveEKE7[1:ekeTimeMean6.shape[1],:],15,cmap=plt.cm.RdBu_r)
#fig.colorbar(cs,ax=ax,ticks=cbarTicks)

ax.contour(x,y,topo7[1:topo6.shape[0]-1,:],15,colors='k')
Q = ax.quiver( X[::plotInterval, ::plotInterval], Y[::plotInterval, ::plotInterval], EP_X_mask7[::plotInterval, ::plotInterval], 
                                          EP_Y_mask7[::plotInterval, ::plotInterval], 
										  pivot='mid', color='k',units='x',scale=0.15e-3,scale_units='x')
#qk = plt.quiverkey(Q, 0.02, 1.1, 3.0*threshold, r'30', coordinates='axes',fontproperties={'weight': 'bold'})
ax.set_xlabel('x (km)')
cbarTicks = np.arange(0.0,61.0,10.0)
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.08, 0.05, 0.82])
fig.colorbar(cs, cax=cbar_ax,ticks=cbarTicks)
cbar_ax.set_title('EKE (kJ.m$^{-3}$)')


fig = plt.figure(4)
ax  = fig.add_subplot(3,1,1)
cs = ax.contourf(x,y,div6,15,cmap=plt.cm.jet)
#ax.contour(x,y,ekeTimeMean6[layerToPlot,1::,1:ekeTimeMean6.shape[2]-1],15,colors='k')
fig.colorbar(cs,ax=ax)

ax  = fig.add_subplot(3,1,2)
cs = ax.contourf(x,y,div5,15,cmap=plt.cm.jet)
#ax.contour(x,y,ekeTimeMean5[layerToPlot,1::,1:ekeTimeMean5.shape[2]-1],15,colors='k')
fig.colorbar(cs,ax=ax)

ax  = fig.add_subplot(3,1,3)
cs = ax.contourf(x,y,div7,15,cmap=plt.cm.jet)
fig.colorbar(cs,ax=ax)


lowerLatitude  = 700.0
upperLatitude = 1350.0


lowerIndex = np.nonzero(y >= lowerLatitude)[0][0]
upperIndex = np.nonzero(y >= upperLatitude)[0][0]

topoTicks = [0.0,1000.0,2000.0]
EKETicks=  [0.0,40.0,80.0,120.0]

#ax.set_yticks([0.2, 0.55, 0.76])
fig = plt.figure(5)
ax1  = fig.add_subplot(5,1,1)

ln1 = ax1.plot(x,densityLayers[layerToPlot]*ekeTimeMean6[layerToPlot,lowerIndex:upperIndex,:].mean(axis=0)*deltaY,'k--',label = 'EKE')
ln2 = ax1.plot(x,densityLayers[layerToPlot]*EPVector_6[0][layerToPlot,lowerIndex:upperIndex,:].mean(axis=0)*deltaY,'k:',label='Wave Flux')
ax1.set_yticks(EKETicks)
ax2 = ax1.twinx()
ln3 = ax2.plot(x,topo6[(lowerIndex+upperIndex)/2,:],'k',label='Topo')
ax2.set_yticks(topoTicks)
# added these three lines
lns = ln1+ln2+ln3
labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, loc='upper left', shadow=True, prop={'size':10})

#ax1.set_xlabel(r"EKE / Wave Flux (J/)")
#ax1.set_ylabel(r"Wave Activity")
#ax2.set_ylabel(r"Height (m)")
ax2.set_xlim(x[0], x[-1])
ax2.set_ylim(0,2500.0)

aax1  = fig.add_subplot(5,1,2)
ln1 = aax1.plot(x,densityLayers[layerToPlot]*ekeTimeMean4[layerToPlot,lowerIndex:upperIndex,:].mean(axis=0)*deltaY,'k--',label = 'EKE')
ln2 = aax1.plot(x,densityLayers[layerToPlot]*EPVector_4[0][layerToPlot,lowerIndex:upperIndex,:].mean(axis=0)*deltaY,'k:',label='Wave Flux')
aax1.set_yticks(EKETicks)

aax2 = aax1.twinx()
ln3 = aax2.plot(x,topo4[(upperIndex),:],'k',label='Topo')
aax2.set_yticks(topoTicks)

aax2.set_xlim(x[0], x[-1])
aax2.set_ylim(0,2500.0)



bx1 = fig.add_subplot(5,1,3)
ln1 = bx1.plot(x,densityLayers[layerToPlot]*ekeTimeMean5[layerToPlot,lowerIndex:upperIndex,:].mean(axis=0)*deltaY,'k--',label = 'EKE')
ln2 = bx1.plot(x,densityLayers[layerToPlot]*EPVector_5[0][layerToPlot,lowerIndex:upperIndex,:].mean(axis=0)*deltaY,'k:',label='Wave Flux')
bx1.set_yticks(EKETicks)

bx2 = bx1.twinx()
ln3 = bx2.plot(x,topo5[(lowerIndex+upperIndex)/2,:],'k',label='Topo')
bx2.set_yticks(topoTicks)

bx2.set_xlim(x[0], x[-1])
bx2.set_ylim(0,2500.0)


bbx1 = fig.add_subplot(5,1,4)
ln1 = bbx1.plot(x,densityLayers[layerToPlot]*ekeTimeMean8[layerToPlot,lowerIndex:upperIndex,:].mean(axis=0)*deltaY,'k--',label = 'EKE')
ln2 = bbx1.plot(x,densityLayers[layerToPlot]*EPVector_8[0][layerToPlot,lowerIndex:upperIndex,:].mean(axis=0)*deltaY,'k:',label='Wave Flux')
bbx1.set_yticks(EKETicks)

bbx2 = bbx1.twinx()
ln3 = bbx2.plot(x,topo8[(lowerIndex+upperIndex)/2,:],'k',label='Topo')
bbx2.set_yticks(topoTicks)

bbx2.set_xlim(x[0], x[-1])
bbx2.set_ylim(0,2500.0)


#bx1.set_xlabel(r"EKE")
#ax1.set_ylabel(r"Wave Activity")
#bx2.set_ylabel(r"Topo")
#bbx2.set_xlim(x[0], x[-1])
#bbx2.set_ylim(0,2500.0)


#ax.plot(x,EPVector_5[0][layerToPlot,lowerIndex:upperIndex,:].mean(axis=0),'b')
#ax.plot(x,ekeTimeMean5[layerToPlot,lowerIndex:upperIndex,:].mean(axis=0),'r')

cx1  = fig.add_subplot(5,1,5)
ln1 = cx1.plot(x,densityLayers[layerToPlot]*ekeTimeMean7[layerToPlot,lowerIndex:upperIndex,:].mean(axis=0)*deltaY,'k--',label = 'EKE')
ln2 = cx1.plot(x,densityLayers[layerToPlot]*EPVector_7[0][layerToPlot,lowerIndex:upperIndex,:].mean(axis=0)*deltaY,'k:',label='Wave Flux')
cx1.set_yticks(EKETicks)

cx2 = cx1.twinx()
ln3 = cx2.plot(x,topo7[(lowerIndex+upperIndex)/2,:],'k',label='Topo')
cx2.set_yticks(topoTicks)

#cx1.set_xlabel(r"EKE")
#ax1.set_ylabel(r"Wave Activity")
#cx2.set_ylabel(r"Topo")
cx2.set_xlim(x[0], x[-1])
cx2.set_ylim(0,2500.0)

import pylab 

pylab.setp(ax1.get_xticklabels(), visible=False)
pylab.setp(bx1.get_xticklabels(), visible=False)
pylab.setp(aax1.get_xticklabels(), visible=False)
pylab.setp(bbx1.get_xticklabels(), visible=False)


bx1.set_ylabel('EKE/Wave Activity Flux (J.m$^{-3}$)', fontsize=14)
bx2.set_ylabel('Topography Height (m)', fontsize=14)

cx1.set_xlabel('x (km)', fontsize=14)



'''
fig = plt.figure(5)
ax  = fig.add_subplot(3,1,1)
cs = ax.contourf(x,y,EPVector_6[0][layerToPlot,...],15,cmap=plt.cm.jet)
ax.contour(x,y,ekeTimeMean6[layerToPlot,1:ekeTimeMean6.shape[1],:],15,colors='k')
fig.colorbar(cs,ax=ax)

ax  = fig.add_subplot(3,1,2)
cs = ax.contourf(x,y,EPVector_5[0][layerToPlot,...],15,cmap=plt.cm.jet)
ax.contour(x,y,ekeTimeMean5[layerToPlot,1:ekeTimeMean6.shape[1],:],15,colors='k')
fig.colorbar(cs,ax=ax)

ax  = fig.add_subplot(3,1,3)
cs = ax.contourf(x,y,EPVector_7[0][layerToPlot,...],15,cmap=plt.cm.jet)
ax.contour(x,y,ekeTimeMean7[layerToPlot,1:ekeTimeMean6.shape[1],:],15,colors='k')
fig.colorbar(cs,ax=ax)
'''


plt.show()

'''
fig = plt.figure(1)
ax  = fig.add_subplot(1,1,1)
ax.contourf(x,y,EPVector_5[2][layerToPlot,...],15,cmap=plt.cm.jet)
ax.contour(x,y,vTimeMean5[layerToPlot,1::,1:vTimeMean5.shape[2]-1],15,colors='k')
ax.contour(x,y,topo5[1:topo5.shape[0]-1,1:topo5.shape[1]-1],15,colors='k')


fig = plt.figure(2)
ax  = fig.add_subplot(1,1,1)
ax.contourf(x,y,EPVector_5[2][layerToPlot,...],15,cmap=plt.cm.jet)
ax.contour(x,y,ekeTimeMean5[layerToPlot,1::,1:ekeTimeMean5.shape[2]-1],15,colors='k')

fig = plt.figure(3)
ax  = fig.add_subplot(1,1,1)
#ax.contourf(x,y,np.sqrt(EPVector_5[0][layerToPlot,...]*EPVector_5[0][layerToPlot,...]),15,cmap=plt.cm.jet)
ax.contourf(x,y,ekeTimeMean5[layerToPlot,1::,1:ekeTimeMean5.shape[2]-1],15,cmap=plt.cm.jet)
ax.contour(x,y,topo5[1:topo5.shape[0]-1,1:topo5.shape[1]-1],15,colors='k')
Q = ax.quiver( X[::22, ::22], Y[::22, ::22], EP_X_mask5[::22, ::22], 
                                          EP_Y_mask5[::22, ::22], 
										  pivot='mid', color='k',units='width',scale=1.0)


fig = plt.figure(4)
ax  = fig.add_subplot(1,1,1)
ax.contourf(x,y,EPVector_6[2][layerToPlot,...],15,cmap=plt.cm.jet)
ax.contour(x,y,vTimeMean6[layerToPlot,1::,1:vTimeMean5.shape[2]-1],15,colors='k')
ax.contour(x,y,topo6[1:topo6.shape[0]-1,1:topo6.shape[1]-1],15,colors='k')


fig = plt.figure(5)
ax  = fig.add_subplot(1,1,1)
ax.contourf(x,y,EPVector_6[2][layerToPlot,...],15,cmap=plt.cm.jet)
ax.contour(x,y,ekeTimeMean6[layerToPlot,1::,1:ekeTimeMean5.shape[2]-1],15,colors='k')


fig = plt.figure(6)
ax  = fig.add_subplot(1,1,1)
#ax.contourf(x,y,np.sqrt(EPVector_5[0][layerToPlot,...]*EPVector_5[0][layerToPlot,...]),15,cmap=plt.cm.jet)
ax.contourf(x,y,ekeTimeMean6[layerToPlot,1::,1:ekeTimeMean5.shape[2]-1],15,cmap=plt.cm.jet)
ax.contour(x,y,topo6[1:topo6.shape[0]-1,1:topo6.shape[1]-1],15,colors='k')
Q = ax.quiver( X[::22, ::22], Y[::22, ::22], EP_X_mask6[::22, ::22], 
                                          EP_Y_mask6[::22, ::22], 
										  pivot='mid', color='k',units='width',scale=1.0)



fig = plt.figure(7)
ax  = fig.add_subplot(1,1,1)
ax.contourf(x,y,EPVector_7[2][layerToPlot,...],15,cmap=plt.cm.jet)
ax.contour(x,y,vTimeMean7[layerToPlot,1::,1:vTimeMean7.shape[2]-1],15,colors='k')
ax.contour(x,y,topo7[1:topo7.shape[0]-1,1:topo7.shape[1]-1],15,colors='k')


fig = plt.figure(8)
ax  = fig.add_subplot(1,1,1)
ax.contourf(x,y,EPVector_7[2][layerToPlot,...],15,cmap=plt.cm.jet)
ax.contour(x,y,ekeTimeMean7[layerToPlot,1::,1:ekeTimeMean7.shape[2]-1],15,colors='k')

fig = plt.figure(9)
ax  = fig.add_subplot(1,1,1)
ax.contourf(x,y,ekeTimeMean7[layerToPlot,1::,1:ekeTimeMean7.shape[2]-1],15,cmap=plt.cm.jet)
ax.contour(x,y,topo7[1:topo7.shape[0]-1,1:topo7.shape[1]-1],15,colors='k')
Q = ax.quiver( X[::22, ::22], Y[::22, ::22], EP_X_mask7[::22, ::22], 
                                          EP_Y_mask7[::22, ::22], 
										  pivot='mid', color='k',units='width',scale=1.0)




fig = plt.figure(10)
ax  = fig.add_subplot(1,1,1)
ax.contourf(x,y,EPVector_8[2][layerToPlot,...],15,cmap=plt.cm.jet)
ax.contour(x,y,vTimeMean8[layerToPlot,1::,1:vTimeMean8.shape[2]-1],15,colors='k')
ax.contour(x,y,topo8[1:topo8.shape[0]-1,1:topo8.shape[1]-1],15,colors='k')


fig = plt.figure(11)
ax  = fig.add_subplot(1,1,1)
ax.contourf(x,y,EPVector_8[2][layerToPlot,...],15,cmap=plt.cm.jet)
ax.contour(x,y,ekeTimeMean8[layerToPlot,1::,1:ekeTimeMean8.shape[2]-1],15,colors='k')

fig = plt.figure(12)
ax  = fig.add_subplot(1,1,1)
ax.contourf(x,y,ekeTimeMean8[layerToPlot,1::,1:ekeTimeMean8.shape[2]-1],15,cmap=plt.cm.jet)
ax.contour(x,y,topo8[1:topo8.shape[0]-1,1:topo7.shape[1]-1],15,colors='k')
Q = ax.quiver( X[::22, ::22], Y[::22, ::22], EP_X_mask8[::22, ::22], 
                                          EP_Y_mask8[::22, ::22], 
										  pivot='mid', color='k',units='width',scale=1.0)


'''







