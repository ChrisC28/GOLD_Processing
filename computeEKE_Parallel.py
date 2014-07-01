import numpy as np
import os
from scipy.io import netcdf as scipy_netcdf
from matplotlib import pyplot as plt
from scipy.signal import firwin, filtfilt
from scipy.fftpack import fft, fftfreq, fftshift, ifft, fft2, ifft2
import multiprocessing
from functools import partial
from time import sleep

def CalcEKERows(iY, nT_total, nT_singleFile, nZ, nX_h, nX_q, fileObjects):

	print iY
	#print nT_total
	#print nT_singleFile
	#print nZ
	#print  nX_h
	#print nX_q
	#return 0,0
	u    = np.zeros((nT_total,nZ,10),dtype='float64')
	v    = np.zeros((nT_total,nZ,10),dtype='float64')
	
	print u.shape
	print v.shape
	#exit()
	u[:,:,:] = 1
	v[:,:,:] = 2
	
	
	#uLowPass  = np.zeros((nT_total,nZ,1,nX_q),dtype='float64')
	#vLowPass  = np.zeros((nT_total,nZ,1,nX_h),dtype='float64')
	#uEddy     = np.zeros((nT_total,nZ,1,nX_q),dtype='float64')
	#vEddy     = np.zeros((nT_total,nZ,1,nX_h),dtype='float64')
	
	#fileCounter  = 0
	#fileObjects = []
	
		
	#for iFile in range(0,len(fileObjects)):

			#Get all the time steps in a single tile 
		#print 'u shape', u[fileCounter*nT_singleFile:(fileCounter+1)*nT_singleFile,:,:].shape
		#print 'v shape', v[fileCounter*nT_singleFile:(fileCounter+1)*nT_singleFile,:,:].shape
		#print 'out shape', uvFile.variables['u'][:,:,iY,nX_q/2:(nX_q/2)+10].shape
		#print 'get some data!'
	#	u[fileCounter*nT_singleFile:(fileCounter+1)*nT_singleFile,:,:]  = fileObjects[iFile].variables['u'][:,:,iY,nX_q/2:(nX_q/2)+10]
	#	v[fileCounter*nT_singleFile:(fileCounter+1)*nT_singleFile,:,:]  = fileObjects[iFile].variables['v'][:,:,iY,nX_h/2:(nX_h/2)+10]
		#print 'got some data!'
		#print u.shape
		#print v.shape
	#	fileCounter = fileCounter + 1
	
	
	
	
	uTimeMean = u.mean(axis=0)
	vTimeMean = v.mean(axis=0)
	print uTimeMean.shape
	uTimeMean = uTimeMean.reshape(nZ*10)
	vTimeMean = vTimeMean.reshape(nZ*10)

	#exit()
	#uTimeMean = uTimeMean.reshape()
	
	return uTimeMean, vTimeMean

def CalcEKEGrid(nT_total, nT_singleFile, nZ_l, nY_h, nY_q, nX_h, nX_q, fileList):
		
	#uTimeMean = np.zeros((nY_h, nX_h),dtype='float64')
	#vTimeMean = np.zeros((nY_h, nX_h),dtype='float64')
	uTimeMean = np.zeros((nZ_l, nY_h, 10),dtype='float64')
	vTimeMean = np.zeros((nZ_l, nY_h, 10),dtype='float64')
	fileObjects = []
	
	print nZ_l
	print nT_total
	for iFileName in fileList[0:2]:
		print iFileName
		fileObjects.append(scipy_netcdf.netcdf_file(iFileName, 'r'))

	 #CalcEKERows(iY, nT_total, nZ_l, nX_h, nX_q, fileList)
	partialCalcRow = partial(CalcEKERows, nT_total=nT_total, nT_singleFile=nT_singleFile, nZ=nZ_l,nX_h=nX_h, nX_q=nX_q,fileObjects=fileObjects)
	
	pool = multiprocessing.Pool(processes = 4)
	print 'starting calculations...'
	
	#parralelOut = pool.map(partialCalcRow, range(0,nY_h))
	parralelOut = pool.map(partialCalcRow, range(0,10))
	print parralelOut
	
	for iFile in range(0,len(fileObjects)):
		fileObjects[iFile].close()	
	sleep(1)
	pool.close()
	pool.join()
	
	uTimeMean = np.asarray(zip(*parralelOut)[0])
	vTimeMean = np.asarray(zip(*parralelOut)[1])
	
	
	print uTimeMean.shape
	print vTimeMean.shape
	uTimeMean = uTimeMean.reshape((nZ_l,10,10))
	vTimeMean = vTimeMean.reshape((nZ_l,10,10))
	

	#uTimeMean, vTimeMean = pool.map(partialCalcRow, range(0,nY_h))
	#uTimeMean = out[0]
	#vTimeMean = out[1]
	return uTimeMean, vTimeMean
	

baseDir = '/home/chris/GOLD/'
fileStartingString = 'ave_prog__'
uvFilesList = []


for file in os.listdir(baseDir):
    if os.path.isfile(file) and file.startswith(fileStartingString):
       uvFilesList.append(baseDir + file)  
uvFilesList = sorted(uvFilesList)

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


print nX_h
print nX_q
print nY_h
print nY_q
deltaT = time_singleFile[1]-time_singleFile[0]
testFile.close()
nT_total  = nT_singleFile * len(uvFilesList)
print nT_singleFile

#print nT_total
#print uvFilesList

tileSizeX = 10
tileSizeY = 10


time = np.zeros(nT_total,dtype='float64')

cutOffPeriod = 90.0
cutOffFreq   = 1.0/cutOffPeriod
#Compute kernal length of the filter
transitionWidth = 2.0/(deltaT)
kernalLength = int(4.0 / transitionWidth)

#Nyquist Frequency
nyquistFrequency = 0.5/deltaT
tapWeights = firwin(kernalLength, cutOffFreq, window='blackman', nyq=nyquistFrequency)
print nT_total
#CalcEKEGrid(nT_total, nZ_l, nY_h, nY_q, nX_h, nX_q, fileList)
uTimeMean, vTimeMean = CalcEKEGrid(nT_total,nT_singleFile, nZ_l, nY_h, nY_q, nX_h, nX_q, uvFilesList)

plt.figure(1)
plt.contourf(uTimeMean[0,:,:],15,cmap=plt.cm.jet)

plt.figure(2)
plt.contourf(yq,vTimeMean[0,:,:],15,cmap=plt.cm.jet)

plt.show()

exit()
'''
tileCounterY = 0
for iTileY in range(0,nY_h,tileSizeY):
	def CalcGrid(y,x):
	outGrid = np.zeros((y.shape[0],x.shape[0]),dtype='float64')
	
	pool =multiprocessing.Pool(processes = 2)
	
	partialCalcRow = partial(calcRow, x=x)

	outGrid = pool.map(partialCalcRow, y)
	pool.close()
	pool.join()

'''
exit()

	
plt.figure(1)
plt.contourf(xq,yh,uTimeMean[0,:,:],15,cmap=plt.cm.jet)

plt.figure(2)
plt.contourf(xh,yq,vTimeMean[0,:,:],15,cmap=plt.cm.jet)

plt.figure(3)
plt.contourf(xq,yh,u2EddyTimeMean[0,:,:],15,cmap=plt.cm.jet)

plt.figure(4)
plt.contourf(xh,yq,v2EddyTimeMean[0,:,:],15,cmap=plt.cm.jet)


	
plt.show()				
exit()	
time = time - time[0]	
deltaT = time[1]-time[0]
u_FFT  = fftshift(fft(u,axis=0))
freqs  = fftshift(fftfreq(nT_total,d=time[1]-time[0]) )	

'''
cutOffPeriod = 90.0
cutOffFreq   = 1.0/cutOffPeriod
#Compute kernal length of the filter
transitionWidth = 2.0/(deltaT)
print transitionWidth
kernalLength = int(4.0 / transitionWidth)

#Nyquist Frequency
print kernalLength
nyquistFrequency = 0.5/deltaT
print nyquistFrequency
tapWeights = firwin(kernalLength, cutOffFreq, window='blackman', nyq=nyquistFrequency)
uLowPass = np.zeros(u.shape,dtype='float64')
vLowPass = np.zeros(v.shape,dtype='float64')

'''
for iZ in range(0,nZ_l):
	uLowPass[:,iZ] = filtfilt(tapWeights,[1.0],u[:,iZ])
	vLowPass[:,iZ] = filtfilt(tapWeights,[1.0],v[:,iZ])

uEddy = u - uLowPass
vEddy = v - vLowPass


uLowPass_FFT  = fftshift(fft(uLowPass,axis=0))


figure1 = plt.figure(1)
ax = figure1.add_subplot(1,1,1)
ax.plot(time,u[:,0],'b')
ax.plot(time,uLowPass[:,0],'r')
ax.plot(time,uEddy[:,0],'g')

figure1 = plt.figure(2)
ax = figure1.add_subplot(1,1,1)
ax.plot(1.0/freqs[nT_total/2::],np.log10(np.abs(u_FFT[nT_total/2::,0]*u_FFT[nT_total/2::,0].conj())),'b')
ax.plot(1.0/freqs[nT_total/2::],np.log10(np.abs(uLowPass_FFT[nT_total/2::,0]*uLowPass_FFT[nT_total/2::,0].conj())),'r')
ax.plot([1.0/cutOffFreq,1.0/cutOffFreq],[-4,3],'k')

plt.show()
