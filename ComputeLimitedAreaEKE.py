import numpy as np
import os
from scipy.io import netcdf as scipy_netcdf
from matplotlib import pyplot as plt
from scipy.signal import firwin, filtfilt, lfilter
from scipy.fftpack import fft, fftfreq, fftshift, ifft, fft2, ifft2


baseDir = '/home/chris/GOLD/IsoTopo4/'
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
deltaX = xh[1]-xh[0]
deltaY = yh[1]-yh[0]
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
 


#Get calculation limits
lowerLatitude  = 700.0
upperLatitude = 2250.0
lowerIndex = np.nonzero(yh >= lowerLatitude)[0][0]
upperIndex = np.nonzero(yh >= upperLatitude)[0][0]

print (upperIndex-lowerIndex) * deltaY

leftLongitude  = 3500.0
rightLongitude =  850.0

leftIndex  = np.nonzero(xh >= leftLongitude)[0][0]
rightIndex = np.nonzero(xh >= rightLongitude)[0][0]
print "left index", leftIndex
print "right index: ", rightIndex

#rightIndex = nX_h-1

lengthY = upperIndex - lowerIndex

if rightIndex > leftIndex:
	lengthX = rightIndex-leftIndex
elif leftIndex > rightIndex:
	lengthX = ((nX_h) - leftIndex) + rightIndex
edgeBuffer = 100.0

time = np.zeros(nT_total,dtype='float64')
u    = np.zeros((nT_total,lengthY+edgeBuffer,lengthX+1+edgeBuffer),dtype='float64')
v    = np.zeros((nT_total,lengthY+1+edgeBuffer,lengthX+edgeBuffer),dtype='float64')


uLowPass  = np.zeros((nT_total,lengthY+edgeBuffer,lengthX+1+edgeBuffer),dtype='float64')
vLowPass  = np.zeros((nT_total,lengthY+1+edgeBuffer,lengthX+edgeBuffer),dtype='float64')

uTransientEddy = np.zeros((nT_total,lengthY+edgeBuffer,lengthX+1+edgeBuffer),dtype='float64')
vTransientEddy = np.zeros((nT_total,lengthY+1+edgeBuffer,lengthX+edgeBuffer),dtype='float64')

EKE = np.zeros((nT_total,lengthY+edgeBuffer,lengthX+edgeBuffer),dtype='float64')

cutOffPeriod = 120.0
cutOffFreq   = 1.0/cutOffPeriod
#Compute kernal length of the filter
transitionWidth = 2.0/(deltaT)
#print transitionWidth
kernalLength = int(4.0 / transitionWidth)

#Nyquist Frequency
nyquistFrequency = 0.5/deltaT
#print nyquistFrequency
tapWeights = firwin(kernalLength, cutOffFreq, window='blackman', nyq=nyquistFrequency)

iLayer = 1
fileCounter  = 0
for iFileName in uvFilesList:
	
	fileObject = scipy_netcdf.netcdf_file(iFileName, 'r')
	time[fileCounter*nT_singleFile:(fileCounter+1)*nT_singleFile] = fileObject.variables['Time'][:]		
	
	
	if rightIndex > leftIndex:		
		u[fileCounter*nT_singleFile:(fileCounter+1)*nT_singleFile,:,:]  = fileObject.variables['u'][:,iLayer,(lowerIndex-edgeBuffer):upperIndex,(leftIndex-edgeBuffer):rightIndex+1]
		v[fileCounter*nT_singleFile:(fileCounter+1)*nT_singleFile,:,:]  = fileObject.variables['v'][:,iLayer,(lowerIndex-edgeBuffer):upperIndex+1,(leftIndex-edgeBuffer):rightIndex]
		
	else:
		print u.shape 
		u[fileCounter*nT_singleFile:(fileCounter+1)*nT_singleFile,:,0:((nX_h)-leftIndex)+edgeBuffer]  = fileObject.variables['u'][:,iLayer,(lowerIndex-edgeBuffer):upperIndex,(leftIndex-edgeBuffer)::]
		u[fileCounter*nT_singleFile:(fileCounter+1)*nT_singleFile,:,((nX_h)-leftIndex)+edgeBuffer::]  = fileObject.variables['u'][:,iLayer,(lowerIndex-edgeBuffer):upperIndex,0:rightIndex+1]
		
		v[fileCounter*nT_singleFile:(fileCounter+1)*nT_singleFile,:,0:((nX_h)-leftIndex)+edgeBuffer]  = fileObject.variables['v'][:,iLayer,(lowerIndex-edgeBuffer):upperIndex+1,(leftIndex-edgeBuffer)::]
		v[fileCounter*nT_singleFile:(fileCounter+1)*nT_singleFile,:,((nX_h)-leftIndex)+edgeBuffer::]  = fileObject.variables['v'][:,iLayer,(lowerIndex-edgeBuffer):upperIndex+1,0:rightIndex]
	
	
	fileCounter = fileCounter + 1
	fileObject.close()
			
#Interpolate to 'h' grid nodes
u_h = 0.5*(u[:,:,1::]+u[:,:,0:-1])
v_h = 0.5*(np.squeeze(v[:,1::,:]+v[:,0:-1,:]))
	
#print "filtering"
uLowPass = lfilter(tapWeights,[1.0],u_h, axis=0)	
vLowPass = lfilter(tapWeights,[1.0],v_h, axis=0)	
#print "finished filtering"

uTransientEddy = u_h - uLowPass
vTransientEddy = v_h - vLowPass
	
	
	
uTransientEddy2 = uTransientEddy*uTransientEddy
vTransientEddy2 = vTransientEddy*vTransientEddy
		
EKE  = 0.5 * (uTransientEddy2 + vTransientEddy2)	
	

postProcessedOutput = scipy_netcdf.netcdf_file('EKE_LimitedArea_IsoTopo4.nc', 'w')

timeDim = postProcessedOutput.createDimension('time', nT_total)
xDim = postProcessedOutput.createDimension('x', lengthX+edgeBuffer)
yDim = postProcessedOutput.createDimension('y', lengthY+edgeBuffer)


timeVar = postProcessedOutput.createVariable('time','f4',('time',))
xVar = postProcessedOutput.createVariable('x','f4',('x',))
yVar = postProcessedOutput.createVariable('y','f4',('y',))

timeVar[:] = time
yVar[:] = yh[lowerIndex-edgeBuffer:upperIndex]
xVar[:] = np.arange(0,(lengthX+edgeBuffer)*deltaX,deltaX)

uVar = postProcessedOutput.createVariable('u','f4',('time','y','x'))
vVar = postProcessedOutput.createVariable('v','f4',('time','y','x'))

uHighPassVar = postProcessedOutput.createVariable('uEddy','f4',('time','y','x'))
vHighPassVar = postProcessedOutput.createVariable('vEddy','f4',('time','y','x'))

EKEVar = postProcessedOutput.createVariable('eke','f4',('time','y','x'))



uVar[:,:,:] = u_h
vVar[:,:,:] = v_h
uHighPassVar[:,:,:] = uTransientEddy
vHighPassVar[:,:,:] = vTransientEddy


EKEVar[:,:,:] = EKE

postProcessedOutput.close()

fig = plt.figure(1)
ax  = fig.add_subplot(1,1,1)
ax.contourf(np.arange(0,(lengthX+edgeBuffer)*deltaX,deltaX),yh[lowerIndex-edgeBuffer:upperIndex],v_h.mean(axis=0),15,cmap=plt.cm.jet)

fig = plt.figure(2)
ax  = fig.add_subplot(1,1,1)
ax.contourf(np.arange(0,(lengthX+edgeBuffer)*deltaX,deltaX),yh[lowerIndex-edgeBuffer:upperIndex],EKE.mean(axis=0),15,cmap=plt.cm.jet)

plt.show()

