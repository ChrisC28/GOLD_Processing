import numpy as np
import os
from scipy.io import netcdf as scipy_netcdf
from matplotlib import pyplot as plt
from scipy.signal import firwin, filtfilt, lfilter
from scipy.fftpack import fft, fftfreq, fftshift, ifft, fft2, ifft2
from netCDF4 import Dataset


baseDir = '/home/chris/GOLD/IsoTopo1/'
outputDirBaseString = 'output'
startYear = 40
endYear   = 70

outputDirs = []

for iYear in range(startYear, endYear):
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
u    = np.zeros((nT_total,nZ_l/2,xqSliceSize+1),dtype='float64')
v    = np.zeros((nT_total,nZ_l/2,2,xhSliceSize),dtype='float64')
e    = np.zeros((nT_total,nZ_i/2+1,xhSliceSize),dtype='float64')


uLowPass  = np.zeros((nT_total,nZ_l/2,xqSliceSize),dtype='float64')
vLowPass  = np.zeros((nT_total,nZ_l/2,xhSliceSize),dtype='float64')

uTransientEddy = np.zeros((nT_total,nZ_l/2,xqSliceSize),dtype='float64')
vTransientEddy = np.zeros((nT_total,nZ_l/2,xhSliceSize),dtype='float64')
eTransientEddy = np.zeros((nT_total,nZ_i/2+1,xqSliceSize),dtype='float64')


uTimeMean = np.zeros((nZ_l/2,nY_h-1,nX_q),dtype='float64')
vTimeMean = np.zeros((nZ_l/2,nY_q-1,nX_h),dtype='float64')
eTimeMean = np.zeros((nZ_i/2+1,nY_q-1,nX_h),dtype='float64')

uZonalMean = np.zeros((nT_total, nZ_l/2,nY_h),dtype='float64')
vZonalMean = np.zeros((nT_total, nZ_l/2,nY_q),dtype='float64')

u2TransientEddyTimeMean = np.zeros((nZ_l/2,nY_h-1,nX_q),dtype='float64')
v2TransientEddyTimeMean = np.zeros((nZ_l/2,nY_h-1,nX_q),dtype='float64')
e2TransientEddyTimeMean = np.zeros((nZ_i/2+1,nY_h-1,nX_q),dtype='float64')

EKE_Transient_TimeMean = np.zeros((nZ_l/2,nY_h-1,nX_q),dtype='float64')




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
		#print '======='
		#print 'iX: ', iX
		#print '======='
		for iFileName in uvFilesList:
			fileObject = scipy_netcdf.netcdf_file(iFileName, 'r')
			#print iFileName
			
			startXIndex = iX*xqSliceSize
			
			if iX < n_Xslice-1:
				#endXIndex = ((iX+1) * xqSliceSize)+1
				xSliceIndicies = range(startXIndex,((iX+1) * xqSliceSize)+1,1)
			else:
				endXIndex = ((iX+1) * xqSliceSize)
				xSliceIndicies = range(startXIndex,((iX+1) * xqSliceSize),1)
				xSliceIndicies.append(0)
			endXIndex = (iX+1)*xqSliceSize
		
			u[fileCounter*nT_singleFile:(fileCounter+1)*nT_singleFile,:,:]  = fileObject.variables['u'][:,::2,iY,xSliceIndicies]
			v[fileCounter*nT_singleFile:(fileCounter+1)*nT_singleFile,:,:,:]  = fileObject.variables['v'][:,::2,iY-1:iY+1,startXIndex:endXIndex]
			e[fileCounter*nT_singleFile:(fileCounter+1)*nT_singleFile,:,:]  = fileObject.variables['e'][:,::2,iY,startXIndex:endXIndex]

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
	
	
	
		uTransientEddy2 = uTransientEddy*uTransientEddy
		vTransientEddy2 = vTransientEddy*vTransientEddy
		eTransientEddy2 = eTransientEddy*eTransientEddy 
		
		#print uTransientEddy2.mean(axis=0).shape
		#print u2TransientEddyTimeMean[:,iY-1,:].shape
#	print np.squeeze(eTransientEddy2.mean(axis=0)).shape
	
		u2TransientEddyTimeMean[:,iY-1,iX*xqSliceSize:(iX+1)*xqSliceSize] = np.squeeze(uTransientEddy2.mean(axis=0))
		v2TransientEddyTimeMean[:,iY-1,iX*xqSliceSize:(iX+1)*xqSliceSize] = np.squeeze(vTransientEddy2.mean(axis=0))
		e2TransientEddyTimeMean[:,iY-1,iX*xqSliceSize:(iX+1)*xqSliceSize] = np.squeeze(eTransientEddy2.mean(axis=0))

	
		EKE_Transient_TimeMean[:,iY-1,iX*xqSliceSize:(iX+1)*xqSliceSize]  = 0.5 * (u2TransientEddyTimeMean[:,iY-1,iX*xqSliceSize:(iX+1)*xqSliceSize] + v2TransientEddyTimeMean[:,iY-1,iX*xqSliceSize:(iX+1)*xqSliceSize])	
	
		
	print "row", iY, 'of ', nY_h, 'done'
	

postProcessedOutput = Dataset('postprocessed_IsoTopo1.nc', 'w', format='NETCDF3_64BIT')

xDim = postProcessedOutput.createDimension('x', size=xh.shape[0])
yDim = postProcessedOutput.createDimension('y', size=yh.shape[0]-1)
zlDim = postProcessedOutput.createDimension('zl', size=zl.shape[0]/2)
ziDim = postProcessedOutput.createDimension('zi', size=zi.shape[0]/2+1)



xVar = postProcessedOutput.createVariable('x','f4',('x',))
yVar = postProcessedOutput.createVariable('y','f4',('y',))
zlVar = postProcessedOutput.createVariable('zl','f4',('zl',))
ziVar = postProcessedOutput.createVariable('zi','f4',('zi',))

xVar[:] = xh
yVar[:] = yh[1::]
zlVar[:] = zl[::2]
ziVar[:] = zi[::2]


uTimeMeanVar = postProcessedOutput.createVariable('u_time_mean','f4',('zl','y','x'))
vTimeMeanVar = postProcessedOutput.createVariable('v_time_mean','f4',('zl','y','x'))
eTimeMeanVar = postProcessedOutput.createVariable('e_time_mean','f4',('zi','y','x'))

EKE_transientVar = postProcessedOutput.createVariable('eke','f4',('zl','y','x'))
e_transientVar = postProcessedOutput.createVariable('e_eddy','f4',('zi','y','x'))

print uTimeMean.shape


uTimeMeanVar[:,:,:] = uTimeMean
vTimeMeanVar[:,:,:] = vTimeMean
eTimeMeanVar[:,:,:] = eTimeMean
e_transientVar[:,:,:] = e2TransientEddyTimeMean
EKE_transientVar[:,:,:] = EKE_Transient_TimeMean

postProcessedOutput.close()


	
plt.figure(1)
plt.contourf(xh,yh[1::],uTimeMean[0,:,:],15,cmap=plt.cm.jet)

plt.figure(2)
plt.contourf(xh,yh[1::],vTimeMean[0,:,:],15,cmap=plt.cm.jet)


plt.figure(5)
plt.contourf(xh,yh[1::],EKE_Transient_TimeMean[0,:,:],15,cmap=plt.cm.jet)

	
plt.show()				
