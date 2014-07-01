import numpy as np
import os
from scipy.io import netcdf as scipy_netcdf
from matplotlib import pyplot as plt
from scipy.signal import firwin, filtfilt, lfilter
from scipy.fftpack import fft, fftfreq, fftshift, ifft, fft2, ifft2
import MontgomeryPot as MP

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
#print densityLayers
#exit()

baseDir = '/home/chris/GOLD/IsoTopo4/'
outputDirBaseString = 'output'
startYear = 39
endYear   = 40

outputDirs = []

for iYear in range(startYear, endYear+1):
	outputDirs.append(outputDirBaseString + "%03d" % iYear + '/')
print outputDirs
aveFileStartingString = 'ave_prog__'
progFileStartingString = 'prog__'
uvFilesList = []
pvFilesList = []
for iYear in range(len(outputDirs)):
	#print baseDir + outputDirs[iYear]
	for file in os.listdir(baseDir + outputDirs[iYear]):
		if os.path.isfile(baseDir + outputDirs[iYear] + file) and file.startswith(aveFileStartingString):
			uvFilesList.append(baseDir + outputDirs[iYear] + file)  

	for file in os.listdir(baseDir + outputDirs[iYear]):
		if os.path.isfile(baseDir + outputDirs[iYear] + file) and file.startswith(progFileStartingString):
			pvFilesList.append(baseDir + outputDirs[iYear] + file) 

uvFilesList = sorted(uvFilesList)
pvFilesList = sorted(pvFilesList)

#print uvFilesList

timeVarName = 'Time'
xhName = 'xh'
xqName = 'xq'
yhName = 'yh'
yqName = 'yq'
zlName = 'zl'
ziName = 'zi'

#get number of time steps

testFile = scipy_netcdf.netcdf_file(uvFilesList[0], 'r')

xh       = testFile.variables[xhName][:]
xq       = testFile.variables[xqName][:]
yh       = testFile.variables[yhName][:]
yq       = testFile.variables[yqName][:]
zl       = testFile.variables[zlName][:]
zi       = testFile.variables[ziName][:]

deltaX  =xq[1]-xq[0]
deltaY  =yq[1]-yq[0]

time_singleFile = testFile.variables[timeVarName][:]
nT_singleFile = time_singleFile.shape[0]
nX_h = xh.shape[0] 
nX_q = xq.shape[0] 
nY_h = yh.shape[0] 
nY_q = yq.shape[0] 
nZ_l = zl.shape[0]
nZ_i = zi.shape[0]

u = testFile.variables['u'][:,:,:,:]
v = testFile.variables['v'][:,:,:,:]
h = testFile.variables['h'][:,:,:,:]
e = testFile.variables['e'][:,:,:,:]

deltaX = xq[1]-xq[0]
deltaY = yq[1]-yq[0]








#dhdy = (h[:,:,2:nY_h,1:nX_h-1] - h[:,:,0:nY_h-2,1:nX_h-1] )/ (2.0*deltaY*1.0e3)
#dhdx = (h[:,:,1:nY_h-1,2:nX_h] - h[:,:,1:nY_h-1,0:nX_h-2] )/ (2.0*deltaX*1.0e3)


#dedy = (e[:,:,2:nY_h,1:nX_h-1] - e[:,:,0:nY_h-2,1:nX_h-1] )/ (2.0*deltaY*1.0e3)
#dedx = (e[:,:,1:nY_h-1,2:nX_h] - e[:,:,1:nY_h-1,0:nX_h-2] )/ (2.0*deltaX*1.0e3)

#dMdx = (montPot[:,:,1:nY_h-1,2:nX_h] - montPot[:,:,1:nY_h-1,0:nX_h-2] )/ (2.0*deltaX*1.0e3)
#dMdy = (montPot[:,:,2:nY_h,1:nX_h-1] - montPot[:,:,0:nY_h-2,1:nX_h-1] )/ (2.0*deltaY*1.0e3)

#vGeo = np.zeros(dMdx.shape,dtype='float64')
#uGeo = np.zeros(dMdy.shape,dtype='float64')

#for iLayer in range(0,nZ_l):
#	vGeo[:,iLayer,:,:] =  1.0/(F0) * dMdx[:,iLayer,:,:]
#	uGeo[:,iLayer,:,:] = -1.0/(F0) * dMdy[:,iLayer,:,:]

uZonalMean = u.mean(axis=-1)
vZonalMean = v.mean(axis=-1)
eZonalMean = e.mean(axis=-1)
hZonalMean = h.mean(axis=-1)


uEddy = np.zeros(u.shape,dtype='float64')
vEddy = np.zeros(v.shape,dtype='float64')
eEddy = np.zeros(e.shape,dtype='float64')
hEddy = np.zeros(h.shape,dtype='float64')

for iT in range(0,nT_singleFile):
	for iZ in range(0,nZ_l):
		for iY in range(0,nY_h):
			uEddy[iT,iZ,iY,:] = u[iT,iZ,iY,:] - uZonalMean[iT,iZ,iY]
			vEddy[iT,iZ,iY,:] = v[iT,iZ,iY,:] - vZonalMean[iT,iZ,iY]
			eEddy[iT,iZ,iY,:] = e[iT,iZ,iY,:] - eZonalMean[iT,iZ,iY]
			hEddy[iT,iZ,iY,:] = h[iT,iZ,iY,:] - hZonalMean[iT,iZ,iY]
			

#np.tile(a,j*i).reshape((i,j,k))
uZonalMean = np.swapaxes(np.tile(uZonalMean,[nX_h]).reshape(nT_singleFile,nZ_l,nX_h,nY_h),2,3)
vZonalMean = np.swapaxes(np.tile(vZonalMean,[nX_h]).reshape(nT_singleFile,nZ_l,nX_h,nY_h),2,3)


montPot = MP.ComputeMontPot(eEddy,hEddy,reducedGravity,densityLayers, RHO0,G0)
#EP_Vector = MP.ComputeEPvector(uEddy[0:10,:,:,:],vEddy[0:10,:,:,:],montPot[0:10,:,:,:], uZonalMean[0:10,:,:,:], vZonalMean[0:10,:,:,:], refThickness, densityLayers, deltaX*1.0e-3, deltaY*1.0e-3, 10.0*G0,F0,RHO0)

#np.save('EP_vector_1.np',EP_Vector[0])
#np.save('EP_vector_2.np',EP_Vector[1])
#np.save('EP_vector_3.np',EP_Vector[2])
#>>> data = np.load('/tmp/123.npz')
EP_1 = np.load('EP_vector_1.npy')
EP_2 = np.load('EP_vector_2.npy')
EP_3 = np.load('EP_vector_3.npy')
print EP_3.shape

deltaT = time_singleFile[1]-time_singleFile[0]
testFile.close()

testFilePV = scipy_netcdf.netcdf_file(pvFilesList[0], 'r')

pv = testFilePV.variables['PV'][:,:,:,:]

testFilePV.close()

layerToPlot = 3
fig = plt.figure(1)
ax = fig.add_subplot(1,1,1)
cs = ax.contourf(xq[1:nX_q-1],yh[1:nY_h-1],EP_3[:,layerToPlot,:,:].mean(axis=0),15,cmap=plt.cm.jet)
fig.colorbar(cs,ax=ax)

fig = plt.figure(2)
ax  = fig.add_subplot(1,1,1)
cs = ax.contourf(xq,yq,(v[:,layerToPlot,:,:]).mean(axis=0),15,cmap=plt.cm.jet)
fig.colorbar(cs,ax=ax)

X,Y = np.meshgrid(xh[1:nX_q-1],yh[1:nY_h-1])

fig = plt.figure(3)
Q = plt.quiver( X[::5, ::5], Y[::5, ::5], EP_1[:,layerToPlot,::5, ::5].mean(axis=0), EP_2[:,layerToPlot,::5, ::5].mean(axis=0), 
pivot='mid', color='k', units='inches' )
qk = plt.quiverkey(Q, 0.5, 0.03, 1, r'$1 \frac{m}{s}$', fontproperties={'weight': 'bold'})
#plt.plot( X[::3, ::3], Y[::3, ::3], 'k.')
#plt.axis([-1, 7, -1, 7])

#l,r,b,t = axis()
#dx, dy = r-l, t-b
#axis([l-0.05*dx, r+0.05*dx, b-0.05*dy, t+0.05*dy])



plt.show()

'''
fig = plt.figure(1)
ax  = fig.add_subplot(1,1,1)
cs = ax.contourf(xq,yq,(pv[:,4,:,:]).mean(axis=0),15,cmap=plt.cm.jet)
fig.colorbar(cs,ax=ax)


fig = plt.figure(2)
ax  = fig.add_subplot(1,1,1)
cs = ax.contourf(xq,yq,(u[:,4,:,:]).mean(axis=0),15,cmap=plt.cm.jet)
fig.colorbar(cs,ax=ax)

fig = plt.figure(3)
ax  = fig.add_subplot(1,1,1)
cs = ax.contourf(xq,yq,(v[:,4,:,:]).mean(axis=0),15,cmap=plt.cm.jet)
fig.colorbar(cs,ax=ax)

plt.show()
'''
'''
meanV = v.mean(axis=0)
pvEddy = pv - pv.mean(axis=0)
PVfluxV = np.zeros(pv.shape,dtype='float64')
for iT in range(0,nT_singleFile):
	
	PVfluxV[iT,:,:,:] = pvEddy[iT,:,:,:] *  meanV
	
fig = plt.figure(1)
ax  = fig.add_subplot(1,1,1)
cs = ax.contourf(xq,yq,(PVfluxV[:,4,:,:]).mean(axis=0),15,cmap=plt.cm.jet)
fig.colorbar(cs,ax=ax)
plt.show()
''' 
