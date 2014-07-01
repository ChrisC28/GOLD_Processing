import numpy as np
import os
from scipy.io import netcdf as scipy_netcdf
from matplotlib import pyplot as plt
from scipy.signal import firwin, filtfilt, lfilter
from scipy.fftpack import fft, fftfreq, fftshift, ifft, fft2, ifft2
import MontgomeryPot as MP

baseDir = '/home/chris/GOLD/IsoTopo7/'
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

uTest=  testFile.variables['u'][:]
print np.mean(uTest[:,4,:,:].mean(axis=-1),axis=0).shape
#plt.figure(1)
#plt.plot(np.mean(uTest[:,4,:,:].mean(axis=-1),axis=0),'b')
#plt.show()



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
e    = np.zeros((nT_singleFile,nZ_i,3,nX_h),dtype='float64')

uEddy = np.zeros((nT_singleFile,nZ_l,3,nX_q),dtype='float64')
vEddy = np.zeros((nT_singleFile,nZ_l,3,nX_h),dtype='float64')
eEddy = np.zeros((nT_singleFile,nZ_i,3,nX_h),dtype='float64')
#hEddy = np.zeros((nT_singleFile,nZ_l,3,nX_h),dtype='float64')


uZonalMean = np.zeros((nT_singleFile,nZ_l,nY_h-2),dtype='float64')
vZonalMean = np.zeros((nT_singleFile,nZ_l,nY_q-2),dtype='float64')
eZonalMean = np.zeros((nT_singleFile,nZ_i,nY_h-2),dtype='float64')

uZonalTimeMean = np.zeros((nZ_l,nY_q-2),dtype='float64')
vZonalTimeMean = np.zeros((nZ_l,nY_q-2),dtype='float64')
eZonalTimeMean = np.zeros((nZ_i,nY_q-2),dtype='float64')


uEddyTimeMean = np.zeros([nZ_l,nY_h-2,nX_q],dtype='float64')
vEddyTimeMean = np.zeros([nZ_l,nY_q-2,nX_h],dtype='float64')
eEddyTimeMean = np.zeros([nZ_i,nY_q-2,nX_h],dtype='float64')

fileCounter  = 0


for iFileName in uvFilesList:
	print iFileName
	fileObject = scipy_netcdf.netcdf_file(iFileName, 'r')

	for iY in range(1,nY_h-1):
		u[:,:,:,:]  = fileObject.variables['u'][:,:,iY-1:iY+2,:]
		v[:,:,:,:]  = fileObject.variables['v'][:,:,iY-1:iY+2,:]
		
		e[:,:,:,:]    = fileObject.variables['e'][:,:,iY-1:iY+2,:]
		
		
		uZonalMean[:,:,iY-1] = u[:,:,1,:].mean(axis=2)
		vZonalMean[:,:,iY-1] = v[:,:,1,:].mean(axis=2)
		eZonalMean[:,:,iY-1] = e[:,:,1,:].mean(axis=2)
		
		
		
		for j in range(0,3):
			
		
			uEddy[:,:,j,:] = u[:,:,j,:] - np.swapaxes(np.tile(u[:,:,j,:].mean(axis=-1),[nX_q]).reshape(nT_singleFile,nX_q,nZ_l),1,2)
			vEddy[:,:,j,:] = v[:,:,j,:] - np.swapaxes(np.tile(v[:,:,j,:].mean(axis=-1),[nX_h]).reshape(nT_singleFile,nX_q,nZ_l),1,2)
			eEddy[:,:,j,:] = e[:,:,j,:] - np.swapaxes(np.tile(e[:,:,j,:].mean(axis=-1),[nX_h]).reshape(nT_singleFile,nX_q,nZ_i),1,2)
			
			#eEddy[:,j,:] = e[:,j,:]   - np.swapaxes(np.tile(e[:,j,:].mean(axis=-1),[nX_h]).reshape(nX_q, nT_singleFile),0,1)
		
		#montPotEddy = MP.ComputeMontPot(eEddy,hEddy,reducedGravity,densityLayers, RHO0,G0)

		
		uEddyTimeMean[:,iY-1,:] = uEddyTimeMean[:,iY-1,:] + uEddy[:,:,1,:].mean(axis=0)
		vEddyTimeMean[:,iY-1,:] = vEddyTimeMean[:,iY-1,:] + vEddy[:,:,1,:].mean(axis=0)
		eEddyTimeMean[:,iY-1,:] = eEddyTimeMean[:,iY-1,:] + eEddy[:,:,1,:].mean(axis=0)
		

		uZonalTimeMean[:,iY-1] = uZonalTimeMean[:,iY-1] + uZonalMean[:,:,iY-1].mean(axis=0)
		vZonalTimeMean[:,iY-1] = vZonalTimeMean[:,iY-1] + vZonalMean[:,:,iY-1].mean(axis=0)
		eZonalTimeMean[:,iY-1] = eZonalTimeMean[:,iY-1] + eZonalMean[:,:,iY-1].mean(axis=0)
		
		
		
		
			
	fileCounter = fileCounter + 1
	
	fileObject.close()	

uFlux = np.zeros([nZ_l,nY_h-2,nX_q],dtype='float64')
vFlux = np.zeros([nZ_l,nY_h-2,nX_q],dtype='float64')


uZTM = np.swapaxes(np.tile(uZonalTimeMean,[nX_q]).reshape(nZ_l,  nX_q, nY_h-2,),1,2)
vZTM = np.swapaxes(np.tile(vZonalTimeMean,[nX_q]).reshape(nZ_l,  nX_q, nY_h-2,),1,2)
eZTM = np.swapaxes(np.tile(eZonalTimeMean,[nX_q]).reshape(nZ_i,  nX_q, nY_h-2,),1,2)
layerToPlot = 3
'''
fig = plt.figure(1)
cs = plt.contourf(eZTM[layerToPlot,...],15,cmap=plt.cm.jet)
fig.colorbar(cs)

fig = plt.figure(2)
cs = plt.contourf(eEddyTimeMean[layerToPlot,...],15,cmap=plt.cm.jet)
fig.colorbar(cs)

fig = plt.figure(3)
cs = plt.contourf(vEddyTimeMean[layerToPlot,...],15,cmap=plt.cm.jet)
fig.colorbar(cs)

fig = plt.figure(4)
cs = plt.contourf((eZTM[1::,...]*vEddyTimeMean)[layerToPlot,:,:],15,cmap=plt.cm.jet)
fig.colorbar(cs)


fig = plt.figure(5)
cs = plt.contourf((vZTM*eEddyTimeMean[1::,...])[layerToPlot,:,:],15,cmap=plt.cm.jet)
fig.colorbar(cs)




plt.show()
'''

uFlux =  (uEddyTimeMean*eEddyTimeMean[1::,...]) # + (uZTM*eEddyTimeMean[1::,...]) + (eZTM[1::,...]*uEddyTimeMean)
vFlux =  (vEddyTimeMean*eEddyTimeMean[1::,...]) # + (vZTM*eEddyTimeMean[1::,...]) + (eZTM[1::,...]*vEddyTimeMean)


uFluxGrad = np.gradient(uFlux,1.0,deltaY*1.0e3,deltaX*1.0e3)
vFluxGrad = np.gradient(vFlux,1.0,deltaY*1.0e3,deltaX*1.0e3)

print uFluxGrad[2][layerToPlot,...].shape
print vFluxGrad[2][layerToPlot,...].shape

#fig = plt.figure(1)
#ax  = fig.add_subplot(1,1,1)
#ax.contourf(xh,yh[1:nY_h-1],-uFluxGrad[2][layerToPlot,...] - vFluxGrad[1][layerToPlot,...],15,cmap=plt.cm.jet)
#plt.show()


#plt.figure(1)
#plt.plot(yq[50:-50],uZonalTimeMean[3,49:-49])
#plt.show()

#for iLayer in range(0,nZ_l):
#	meanEKE[iLayer,:,:]  = densityLayers[iLayer]* refThickness[iLayer]*meanEKE[iLayer,:,:] 
#	meanPE[iLayer,:,:] = densityLayers[iLayer] * reducedGravity[iLayer] * meanPE[iLayer,:,:]

postProcessedOutput = scipy_netcdf.netcdf_file('standing_fluxes_isoTopo7.nc', 'w')

xDim = postProcessedOutput.createDimension('x', xh.shape[0])
yDim = postProcessedOutput.createDimension('y', yh.shape[0]-2)
zlDim = postProcessedOutput.createDimension('zl', zl.shape[0])
ziDim = postProcessedOutput.createDimension('zi', zi.shape[0])



xVar = postProcessedOutput.createVariable('x','f4',('x',))
yVar = postProcessedOutput.createVariable('y','f4',('y',))
zlVar = postProcessedOutput.createVariable('zl','f4',('zl',))
ziVar = postProcessedOutput.createVariable('zi','f4',('zi',))

xVar[:] = xh
yVar[:] = yh[1:nY_h-1]
zlVar[:] = zl
ziVar[:] = zi


uZonalMeanVar = postProcessedOutput.createVariable('uMean','f4',('zl','y'))
vZonalMeanVar = postProcessedOutput.createVariable('vMean','f4',('zl','y'))
eZonalMeanVar = postProcessedOutput.createVariable('eMean','f4',('zi','y'))

uEddyTimeMeanVar = postProcessedOutput.createVariable('uEddyTimeMean','f4',('zl','y','x'))
vEddyTimeMeanVar = postProcessedOutput.createVariable('vEddyTimeMean','f4',('zl','y','x'))
eEddyTimeMeanVar = postProcessedOutput.createVariable('eEddyTimeMean','f4',('zi','y','x'))

statEddyFluxUVar = postProcessedOutput.createVariable('statEddyFluxU','f4',('zl','y','x'))
statEddyFluxVVar = postProcessedOutput.createVariable('statEddyFluxV','f4',('zl','y','x'))

uZonalMeanVar[:,:] = uZonalTimeMean
vZonalMeanVar[:,:] = vZonalTimeMean
eZonalMeanVar[:,:] = eZonalTimeMean

uEddyTimeMeanVar[:,:,:] = uEddyTimeMean
vEddyTimeMeanVar[:,:,:] = vEddyTimeMean
eEddyTimeMeanVar[:,:,:] = eEddyTimeMean

statEddyFluxUVar[:,:,:] = uFlux
statEddyFluxVVar[:,:,:] = vFlux





postProcessedOutput.close()








