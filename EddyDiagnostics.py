import numpy as np
import os
from scipy.io import netcdf as scipy_netcdf
from matplotlib import pyplot as plt
from scipy.signal import firwin, filtfilt, lfilter
from scipy.fftpack import fft, fftfreq, fftshift, ifft, fft2, ifft2
import MontgomeryPot as MP

baseDir = '/home/chris/GOLD/RandTopo1/'
outputDirBaseString = 'output'
startYear = 39
endYear   = 59

outputDirs = []
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


for iYear in range(startYear, endYear+1):
	outputDirs.append(outputDirBaseString + "%03d" % iYear + '/')
print outputDirs
fileStartingString = 'ave_prog__'
uvFilesList = []

for iYear in range(len(outputDirs)):
	#print baseDir + outputDirs[iYear]
	for file in os.listdir(baseDir + outputDirs[iYear]):
		if os.path.isfile(baseDir + outputDirs[iYear] + file) and file.startswith(fileStartingString):
			uvFilesList.append(baseDir + outputDirs[iYear] + file)  
uvFilesList = sorted(uvFilesList)
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

time_singleFile = testFile.variables[timeVarName][:]
nT_singleFile = time_singleFile.shape[0]
nX_h = xh.shape[0] 
nX_q = xq.shape[0] 
nY_h = yh.shape[0] 
nY_q = yq.shape[0] 
nZ_l = zl.shape[0]
nZ_i = zi.shape[0]

deltaX = xq[1]-xq[0]
deltaY = yq[1]-yq[0]

deltaT = time_singleFile[1]-time_singleFile[0]
testFile.close()
nT_total  = nT_singleFile * len(uvFilesList)

n_Xslice   = 1
xqSliceSize = nX_q / n_Xslice
xhSliceSize = nX_h / n_Xslice
 



time = np.zeros(nT_total,dtype='float64')

u    = np.zeros((nT_singleFile,nZ_l,3,nX_q),dtype='float64')
v    = np.zeros((nT_singleFile,nZ_l,3,nX_h),dtype='float64')
e    = np.zeros((nT_singleFile,3,nX_h),dtype='float64')
h    = np.zeros((nT_singleFile,nZ_l,3,nX_h),dtype='float64')

uEddy = np.zeros((nT_singleFile,nZ_l,3,nX_q),dtype='float64')
vEddy = np.zeros((nT_singleFile,nZ_l,3,nX_h),dtype='float64')
eEddy = np.zeros((nT_singleFile,3,nX_h),dtype='float64')
hEddy = np.zeros((nT_singleFile,nZ_l,3,nX_h),dtype='float64')


uZonalMean = np.zeros((nT_singleFile,nZ_l,nY_h-2),dtype='float64')
vZonalMean = np.zeros((nT_singleFile,nZ_l,nY_q-2),dtype='float64')
eZonalMean = np.zeros((nT_singleFile,nZ_i,nY_h-2),dtype='float64')
hZonalMean = np.zeros((nT_singleFile,nZ_l,nY_h-2),dtype='float64')


uEddyTimeMean = np.zeros([nZ_l,nY_h-2,nX_q],dtype='float64')
vEddyTimeMean = np.zeros([nZ_l,nY_q-2,nX_h],dtype='float64')
#mEddyTimeMean = np.zeros([nZ_l,nY_q-2,nX_h],dtype='float64')
meanEKE = np.zeros([nZ_l,nY_h-2,nX_h],dtype='float64')
meanPE  = np.zeros([nZ_l,nY_h-2,nX_h],dtype='float64')

EPvectorTimeMean_X = np.zeros([nZ_l,nY_h-2,nX_h-2])
EPvectorTimeMean_Y = np.zeros([nZ_l,nY_h-2,nX_h-2])
EPvectorTimeMean_Z = np.zeros([nZ_l,nY_h-2,nX_h-2])


fileCounter  = 0


for iFileName in uvFilesList:
	print iFileName
	fileObject = scipy_netcdf.netcdf_file(iFileName, 'r')

	for iY in range(1,nY_h-1):
		#print iY
		u[:,:,:,:]  = fileObject.variables['u'][:,:,iY-1:iY+2,:]
		v[:,:,:,:]  = fileObject.variables['v'][:,:,iY-1:iY+2,:]
		h[:,:,:,:]  = fileObject.variables['h'][:,:,iY-1:iY+2,:]
		
		#Get upper layer interface height for the computation of the 
		#Montgomery Potentials
		
		e[:,:,:]    = fileObject.variables['e'][:,0,iY-1:iY+2,:]
		for iLayer in range(0,nZ_l):
			h[:,iLayer,:,:] = h[:,iLayer,:] - refThickness[iLayer]
		
		
		uZonalMean[:,:,iY-1] = u[:,:,1,:].mean(axis=-1)
		vZonalMean[:,:,iY-1] = v[:,:,1,:].mean(axis=-1)
		hZonalMean[:,:,iY-1] = h[:,:,1,:].mean(axis=-1)
		
		for j in range(0,3):
			
		
			uEddy[:,:,j,:] = u[:,:,j,:] - np.swapaxes(np.tile(u[:,:,j,:].mean(axis=-1),[nX_q]).reshape(nT_singleFile,nX_q,nZ_l),1,2)
			vEddy[:,:,j,:] = v[:,:,j,:] - np.swapaxes(np.tile(v[:,:,j,:].mean(axis=-1),[nX_h]).reshape(nT_singleFile,nX_q,nZ_l),1,2)
			hEddy[:,:,j,:] = h[:,:,j,:] - np.swapaxes(np.tile(h[:,:,j,:].mean(axis=-1),[nX_h]).reshape(nT_singleFile,nX_q,nZ_l),1,2)
			
			eEddy[:,j,:] = e[:,j,:]   - np.swapaxes(np.tile(e[:,j,:].mean(axis=-1),[nX_h]).reshape(nX_q, nT_singleFile),0,1)
		
		montPotEddy = MP.ComputeMontPot(eEddy,hEddy,reducedGravity,densityLayers, RHO0,G0)

		
		uEddyTimeMean[:,iY-1,:] = uEddyTimeMean[:,iY-1,:] + uEddy[:,:,1,:].mean(axis=0)
		vEddyTimeMean[:,iY-1,:] = vEddyTimeMean[:,iY-1,:] + vEddy[:,:,1,:].mean(axis=0)
		#mEddyTimeMean[:,iY-1,:] = mEddyTimeMean[:,iY-1,:] + montPotEddy[:,:,1,:].mean(axis=0)
		
		EP_Vector = MP.ComputeEPvector(uEddy,vEddy,montPotEddy, np.swapaxes(np.tile(uZonalMean[:,:,iY-1],[nX_q]).reshape(nT_singleFile,nX_q,nZ_l),1,2), 
										np.swapaxes(np.tile(vZonalMean[:,:,iY-1],[nX_h]).reshape(nT_singleFile,nX_q,nZ_l),1,2), 
										refThickness, densityLayers, deltaX*1.0e3, deltaY*1.0e3, 10.0*G0,F0,RHO0)
		
		EPvectorTimeMean_X[:,iY-1,:] = EPvectorTimeMean_X[:,iY-1,:] + EP_Vector[0].mean(axis=0)
		EPvectorTimeMean_Y[:,iY-1,:] = EPvectorTimeMean_Y[:,iY-1,:] + EP_Vector[1].mean(axis=0)
		EPvectorTimeMean_Z[:,iY-1,:] = EPvectorTimeMean_Z[:,iY-1,:] + EP_Vector[2].mean(axis=0)
		
		
		
		meanEKE[:,iY-1,:] = meanEKE[:,iY-1,:] + (0.5*( uEddy[:,:,1,:]*uEddy[:,:,1,:] + vEddy[:,:,1,:]*vEddy[:,:,1,:] )).mean(axis=0)
		meanPE[:,iY-1,:]  = meanPE[:,iY-1,:]  + (0.5*( hEddy[:,:,1,:]*hEddy[:,:,1,:] )).mean(axis=0)	
		
		
		
		
			
	fileCounter = fileCounter + 1
	
	fileObject.close()	


for iLayer in range(0,nZ_l):
	meanEKE[iLayer,:,:]  = densityLayers[iLayer]* refThickness[iLayer]*meanEKE[iLayer,:,:] 
	meanPE[iLayer,:,:] = densityLayers[iLayer] * reducedGravity[iLayer] * meanPE[iLayer,:,:]

postProcessedOutput = scipy_netcdf.netcdf_file('waveActivity_Zonal_randTopoSmooth1.nc', 'w')

xDim = postProcessedOutput.createDimension('x', xh.shape[0]-2)
yDim = postProcessedOutput.createDimension('y', yh.shape[0]-2)
zlDim = postProcessedOutput.createDimension('zl', zl.shape[0])



xVar = postProcessedOutput.createVariable('x','f4',('x',))
yVar = postProcessedOutput.createVariable('y','f4',('y',))
zlVar = postProcessedOutput.createVariable('zl','f4',('zl',))

xVar[:] = xh[1:nX_h-1]
yVar[:] = yh[1:nY_h-1]
zlVar[:] = zl


EPXVar = postProcessedOutput.createVariable('EP_X','f4',('zl','y','x'))
EPYVar = postProcessedOutput.createVariable('EP_Y','f4',('zl','y','x'))
EPZVar = postProcessedOutput.createVariable('EP_Z','f4',('zl','y','x'))

PE_Var = postProcessedOutput.createVariable('PE','f4',('zl','y','x'))



EPXVar[:,:,:] = EPvectorTimeMean_X
EPYVar[:,:,:] = EPvectorTimeMean_Y
EPZVar[:,:,:] = EPvectorTimeMean_Z
PE_Var[:,:,:] = meanPE[:,:,1:nX_h-1]

postProcessedOutput.close()







layerToPlot = 3

normEP = np.sqrt(EPvectorTimeMean_X**2+EPvectorTimeMean_Y**2)

maskArray = normEP[layerToPlot,:,:] < 1.0*normEP[layerToPlot,:,:].std()
EP_X_mask = np.ma.masked_array(EPvectorTimeMean_X[layerToPlot,:,:], mask=maskArray)
EP_Y_mask = np.ma.masked_array(EPvectorTimeMean_Y[layerToPlot,:,:], mask=maskArray)


fig = plt.figure(1)
ax = fig.add_subplot(1,1,1)
cs = ax.contourf(xq,yh[1:nY_h-1],uEddyTimeMean[layerToPlot,:,:],15,cmap=plt.cm.jet)
fig.colorbar(cs,ax=ax)


fig = plt.figure(2)
ax = fig.add_subplot(1,1,1)
cs = ax.contourf(xh,yq[1:nY_q-1],vEddyTimeMean[layerToPlot,:,:],15,cmap=plt.cm.jet)
fig.colorbar(cs,ax=ax)


fig = plt.figure(3)
ax = fig.add_subplot(1,1,1)
cs = ax.contourf(xh[1:nX_h-1],yq[1:nY_q-1],EPvectorTimeMean_Z[layerToPlot,:,:],15,cmap=plt.cm.jet)
fig.colorbar(cs,ax=ax)

X,Y = np.meshgrid(xh,yh)

#fig = plt.figure(4)
Q = ax.quiver( X[::20, ::20], Y[::20, ::20], EP_X_mask[::20, ::20], 
                                          EP_Y_mask[::20, ::20], 
										  pivot='mid', color='k',units='width',scale=1.0)
#ax.set_xlim([xh[0],xh[-1]])
#ax.set_xlim([yh[0],yh[-1]])

#qk = plt.quiverkey(Q, 0.5, 0.03, 1, r'$1 \frac{m}{s}$', fontproperties={'weight': 'bold'})

fig = plt.figure(5)
ax = fig.add_subplot(1,1,1)
cs = ax.contourf(xh[1:nX_h-1],yq[1:nY_q-1], EPvectorTimeMean_X[layerToPlot,:,:],15,cmap=plt.cm.jet)
fig.colorbar(cs,ax=ax)
fig = plt.figure(6)
ax = fig.add_subplot(1,1,1)
cs = ax.contourf(xh[1:nX_h-1],yq[1:nY_q-1], EPvectorTimeMean_Y[layerToPlot,:,:],15,cmap=plt.cm.jet)
fig.colorbar(cs,ax=ax)


plt.show()

