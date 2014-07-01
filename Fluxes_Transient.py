import numpy as np
import os
from scipy.io import netcdf as scipy_netcdf
from matplotlib import pyplot as plt
from scipy.signal import firwin, filtfilt, lfilter
from scipy.fftpack import fft, fftfreq, fftshift, ifft, fft2, ifft2


baseDir = '/home/chris/GOLD/RandTopoSmooth/'
outputDirBaseString = 'output'
startYear = 39
endYear   = 59

outputDirs = []

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


deltaT = time_singleFile[1]-time_singleFile[0]
testFile.close()
nT_total  = nT_singleFile * len(uvFilesList)

print nT_singleFile
print nT_total
#print uvFilesList
tileSizeX = 10
tileSizeY = 10

n_Xslice   = 32
xqSliceSize = nX_q / n_Xslice
xhSliceSize = nX_h / n_Xslice
 


time = np.zeros(nT_total,dtype='float64')
u    = np.zeros((nT_total,nZ_l,xqSliceSize+1),dtype='float64')
v    = np.zeros((nT_total,nZ_l,2,xhSliceSize),dtype='float64')
e    = np.zeros((nT_total,nZ_i,xhSliceSize),dtype='float64')


uLowPass  = np.zeros((nT_total,nZ_l,xqSliceSize),dtype='float64')
vLowPass  = np.zeros((nT_total,nZ_l,xhSliceSize),dtype='float64')

uTransientEddy = np.zeros((nT_total,nZ_l,xqSliceSize),dtype='float64')
vTransientEddy = np.zeros((nT_total,nZ_l,xhSliceSize),dtype='float64')
eTransientEddy = np.zeros((nT_total,nZ_i,xqSliceSize),dtype='float64')


uTimeMean = np.zeros((nZ_l,nY_h-1,nX_q),dtype='float64')
vTimeMean = np.zeros((nZ_l,nY_q-1,nX_h),dtype='float64')
eTimeMean = np.zeros((nZ_i,nY_q-1,nX_h),dtype='float64')


uTotalFluxTimeMean = np.zeros((nZ_l,nY_h-1,nX_q),dtype='float64')
vTotalFluxTimeMean = np.zeros((nZ_l,nY_h-1,nX_q),dtype='float64')

uEddyFluxTimeMean = np.zeros((nZ_l,nY_h-1,nX_q),dtype='float64')
vEddyFluxTimeMean = np.zeros((nZ_l,nY_h-1,nX_q),dtype='float64')
eEddyTimeMean     = np.zeros((nZ_i,nY_h-1,nX_q),dtype='float64')

cutOffPeriod = 180.0
cutOffFreq   = 1.0/cutOffPeriod
#Compute kernal length of the filter
transitionWidth = 2.0/(deltaT)
#print transitionWidth
kernalLength = int(4.0 / transitionWidth)

#Nyquist Frequency
nyquistFrequency = 0.5/deltaT
#print nyquistFrequency
tapWeights = firwin(kernalLength, cutOffFreq, window='blackman', nyq=nyquistFrequency)


for iY in range(1,nY_h):
		
	for iX in range(0,n_Xslice):

		fileCounter  = 0
		for iFileName in uvFilesList:
			
			fileObject = scipy_netcdf.netcdf_file(iFileName, 'r')
			
			startXIndex = iX*xqSliceSize
			
			if iX < n_Xslice-1:
				#endXIndex = ((iX+1) * xqSliceSize)+1
				xSliceIndicies = range(startXIndex,((iX+1) * xqSliceSize)+1,1)
			else:
				endXIndex = ((iX+1) * xqSliceSize)
				xSliceIndicies = range(startXIndex,((iX+1) * xqSliceSize),1)
				xSliceIndicies.append(0)
			endXIndex = (iX+1)*xqSliceSize
		
			u[fileCounter*nT_singleFile:(fileCounter+1)*nT_singleFile,:,:]  = fileObject.variables['u'][:,:,iY,xSliceIndicies]
			v[fileCounter*nT_singleFile:(fileCounter+1)*nT_singleFile,:,:,:]  = fileObject.variables['v'][:,:,iY-1:iY+1,startXIndex:endXIndex]
			e[fileCounter*nT_singleFile:(fileCounter+1)*nT_singleFile,:,:]  = fileObject.variables['e'][:,:,iY,startXIndex:endXIndex]

			fileCounter = fileCounter + 1
			fileObject.close()
			
		#Interpolate to 'h' grid nodes
		u_h = 0.5*(u[:,:,1::]+u[:,:,0:-1])
		v_h = 0.5*(np.squeeze(v[:,:,1::,:]+v[:,:,0:-1,:]))
	
		#print "filtering"
		uLowPass = lfilter(tapWeights,[1.0],u_h, axis=0)	
		vLowPass = lfilter(tapWeights,[1.0],v_h, axis=0)	
		eLowPass = lfilter(tapWeights,[1.0],e, axis=0)	
		#print "finished filtering"
		
		uTransientEddy = u_h - uLowPass
		vTransientEddy = v_h - vLowPass
		eTransientEddy = e - eLowPass
		
		#print vTimeMean[:,iY-1,iX*xqSliceSize:(iX+1)*xqSliceSize].shape
		#print np.squeeze(v_h.mean(axis=0)).shape
		uTimeMean[:,iY-1,iX*xqSliceSize:(iX+1)*xqSliceSize] = np.squeeze(u_h.mean(axis=0))
		vTimeMean[:,iY-1,iX*xqSliceSize:(iX+1)*xqSliceSize] = np.squeeze(v_h.mean(axis=0))
		eTimeMean[:,iY-1,iX*xqSliceSize:(iX+1)*xqSliceSize] = np.squeeze(e.mean(axis=0))
		
		#print uTransientEddy.shape
		#print eTransientEddy[1::,...].shape
	
		uFlux = uTransientEddy*eTransientEddy[:,1::,...]
		vFlux = vTransientEddy*eTransientEddy[:,1::,...]
		
		uEddyFluxTimeMean[:,iY-1,iX*xqSliceSize:(iX+1)*xqSliceSize] = np.squeeze(uFlux.mean(axis=0))
		vEddyFluxTimeMean[:,iY-1,iX*xqSliceSize:(iX+1)*xqSliceSize] = np.squeeze(vFlux.mean(axis=0))
		
		uFlux = u_h * e[:,1::,...]
		vFlux = v_h * e[:,1::,...]
		
		uTotalFluxTimeMean[:,iY-1,iX*xqSliceSize:(iX+1)*xqSliceSize] = np.squeeze(uFlux.mean(axis=0))
		vTotalFluxTimeMean[:,iY-1,iX*xqSliceSize:(iX+1)*xqSliceSize] = np.squeeze(vFlux.mean(axis=0))
		eEddyTimeMean[:,iY-1,iX*xqSliceSize:(iX+1)*xqSliceSize] = np.squeeze(eTransientEddy.mean(axis=0))
	print "row", iY, 'of ', nY_h, 'done'
	

postProcessedOutput = scipy_netcdf.netcdf_file('transient_fluxes_randTopoSmooth.nc', 'w')

xDim = postProcessedOutput.createDimension('x', xh.shape[0])
yDim = postProcessedOutput.createDimension('y', yh.shape[0]-1)
zlDim = postProcessedOutput.createDimension('zl', zl.shape[0])
ziDim = postProcessedOutput.createDimension('zi', zi.shape[0])



xVar = postProcessedOutput.createVariable('x','f4',('x',))
yVar = postProcessedOutput.createVariable('y','f4',('y',))
zlVar = postProcessedOutput.createVariable('zl','f4',('zl',))
ziVar = postProcessedOutput.createVariable('zi','f4',('zi',))

xVar[:] = xh
yVar[:] = yh[1::]
zlVar[:] = zl
ziVar[:] = zi


uTimeMeanVar = postProcessedOutput.createVariable('u_time_mean','f4',('zl','y','x'))
vTimeMeanVar = postProcessedOutput.createVariable('v_time_mean','f4',('zl','y','x'))
eTimeMeanVar = postProcessedOutput.createVariable('e_time_mean','f4',('zi','y','x'))
eEddyTimeMeanVar = postProcessedOutput.createVariable('e_eddy_time_mean','f4',('zi','y','x'))


uFluxTotalVar = postProcessedOutput.createVariable('uTotalFlux','f4',('zl','y','x'))
vFluxTotalVar = postProcessedOutput.createVariable('vTotalFlux','f4',('zl','y','x'))





uFluxEddyVar = postProcessedOutput.createVariable('uEddyFlux','f4',('zl','y','x'))
vFluxEddyVar = postProcessedOutput.createVariable('vEddyFlux','f4',('zl','y','x'))

print uTimeMean.shape


uTimeMeanVar[:,:,:] = uTimeMean
vTimeMeanVar[:,:,:] = vTimeMean
eTimeMeanVar[:,:,:] = eTimeMean
eEddyTimeMeanVar[:,:,:] = eEddyTimeMean
uFluxTotalVar[:,:,:] = uTotalFluxTimeMean
vFluxTotalVar[:,:,:] = vTotalFluxTimeMean


uFluxEddyVar[:,:,:] = uEddyFluxTimeMean
vFluxEddyVar[:,:,:] = vEddyFluxTimeMean

postProcessedOutput.close()


layerToPlot = 4
plt.figure(1)
plt.contourf(xh,yh[1::],uTimeMean[layerToPlot,:,:],15,cmap=plt.cm.jet)

plt.figure(2)
plt.contourf(xh,yh[1::],vTimeMean[layerToPlot,:,:],15,cmap=plt.cm.jet)


plt.figure(3)
plt.contourf(xh,yh[1::],uFluxVar[layerToPlot,:,:],15,cmap=plt.cm.jet)

plt.figure(4)
plt.contourf(xh,yh[1::],vFluxVar[layerToPlot,:,:],15,cmap=plt.cm.jet)
plt.show()				
