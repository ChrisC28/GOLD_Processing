import numpy as np
import os
from scipy.io import netcdf as scipy_netcdf
from matplotlib import pyplot as plt
from scipy.signal import firwin, filtfilt, lfilter
from scipy.fftpack import fft, fftfreq, fftshift, ifft, fft2, ifft2, fftn, ifftn
from scipy.linalg import lstsq
from scipy.ndimage.filters import gaussian_filter
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
		
		
	
def FitPlane(inputData,x,y):
	
	
	inputRearranged  = np.reshape(inputData, (x.shape[0]*y.shape[0]))
	X_grid, Y_grid = np.meshgrid(x, y)	

	X_grid = np.reshape(X_grid,x.shape[0]*y.shape[0])
	Y_grid = np.reshape(Y_grid,x.shape[0]*y.shape[0])

	A = np.column_stack([X_grid, Y_grid, np.ones_like(X_grid)])
	planeParameters = lstsq(A, inputRearranged, cond=None, overwrite_a=False, overwrite_b=False)
	print planeParameters[0]

	X_grid = np.reshape(X_grid,(y.shape[0],x.shape[0]) )
	Y_grid = np.reshape(Y_grid,(y.shape[0],x.shape[0]) )

	#topoGrid  = np.reshape(topoGrid, (nYgrid,nXgrid))
	meanPlane = (planeParameters[0][0]*X_grid) + (planeParameters[0][1]*Y_grid) + planeParameters[0][2]
	return meanPlane

def DesignGaussFilter(cutoffWavenumber,valueAtCutoff, waveNumbers):
	#print cutoffWavenumber
	
	#length = np.nonzero(waveNumbers[waveNumbers.shape[0]/2::] >= cutoffWavenumber)[0][0]
	
	
	sigma_k =  np.sqrt( - (0.5 *cutoffWavenumber*cutoffWavenumber)  / np.log(valueAtCutoff))
	return sigma_k
	
	
	
	
	

timeVarName = 'time'
xName = 'x'
yName = 'y'

#get number of time steps
fileName = "/home/chris/GOLD/IsoTopo_4_Processing/EKE_LimitedArea_IsoTopo4.nc"
dataFile = scipy_netcdf.netcdf_file(fileName, 'r')

time   = dataFile.variables[timeVarName][:]
x      = dataFile.variables[xName][:]
y      = dataFile.variables[yName][:]

EKE    = dataFile.variables['eke'][:]


deltaX  = x[1]-x[0]
deltaY  = y[1]-y[0]
nX = x.shape[0]
nY = y.shape[0]

print y[-1]-y[0]
nFilters = 7


lowerLatitude = 800.0
upperLatitude = 1650.0
lowerIndex = np.nonzero(y >= lowerLatitude)[0][0]
upperIndex = np.nonzero(y >= upperLatitude)[0][0]

nyquistFrequencyX = 0.5/deltaX
nyquistFrequencyY = 0.5/deltaY
kernalLength = 40


xWaveNumbers = fftshift(fftfreq(x.shape[0],d=deltaX))
yWaveNumbers = fftshift(fftfreq(y.shape[0],d=deltaY))

EKE_Plane = FitPlane(EKE[0,:,:],x,y)

'''
plt.figure(1)
plt.contourf(x,y,EKE[0,:,:],15,cmap=plt.cm.jet)
plt.colorbar()

plt.figure(2)
plt.contourf(x,y,EKE_Plane,15,cmap=plt.cm.jet)
plt.colorbar()
plt.figure(3)
plt.contourf(x,y,EKE[0,:,:]-EKE_Plane,15,cmap=plt.cm.jet)
plt.colorbar()
plt.show()
exit()
'''

cutoffsX = []
cutoffsY = []
EKE_IntegratedX = []
EKE_IntegratedY = []

EKE_NEW_X = np.zeros(EKE[0,:,:].shape,dtype='float64')
EKE_NEW_Y = np.zeros(EKE.shape,dtype='float64')


bufferLength = 2*kernalLength

cutoffs = []
cutoffs = []
filterBandWith = 3
upperCutoffY = nY/2

print nX
for iBand in range(0,int(yWaveNumbers[nY/2::].shape[0]/filterBandWith)):
	upperCutoffY = upperCutoffY+filterBandWith
#	print upperCutoffX
	
	cutoffs.append(yWaveNumbers[upperCutoffY-1])
	#print 1.0/cutoffs[-1]

kGrid,lGrid = np.meshgrid(xWaveNumbers, yWaveNumbers)


cutoffs.reverse()
EKE = EKE.astype('float64')

EKE_Int_FromFFT = []

EKE_NEW_X[:,:] = EKE[0,:,:] -EKE_Plane

#fft_EKE = fftshift(fftn(EKE,axes=[2,1]))

#fft_EKE = fftshift(fftn((EKE[0,:,:]-EKE_Plane) - np.mean(EKE[0,:,:]-EKE_Plane),axes=[1,0]))
'''
valueAtCutoff = 0.01
for iFilter in range(10,len(cutoffs)):
	print 1.0/cutoffs[iFilter]
	#EKE_BP_FFT = np.zeros(EKE[0,:,:].shape,dtype='float64')
	
	
	sigma = DesignGaussFilter(cutoffs[iFilter],valueAtCutoff,yWaveNumbers)
	
	
	#print sigma
	filterKernal= fftshift(np.exp( -((kGrid*kGrid) + (lGrid*lGrid))/(2.0*sigma*sigma)))
	EKE_LP = ifftn((fftn(EKE_NEW,axes=[1,0])*filterKernal),axes=[1,0])
	
	EKE_HP = EKE_NEW-EKE_LP
	
	
	#plt.figure(1)
	#plt.plot(xWaveNumbers,filterKernal[nY/2,:])
	#
	#plt.show()
	
	
	
	#print sigma
	#EKE_LP = gaussian_filter(EKE_NEW, sigma, mode='mirror')
	#EKE_HP = EKE_NEW-EKE_LP
	
	
	plt.figure(1)
	plt.contourf(x,y,EKE_NEW,15,cmap=plt.cm.jet)
	plt.colorbar()
	
	plt.figure(2)
	plt.contourf(x,y,EKE_LP,15,cmap=plt.cm.jet)
	plt.colorbar()
	
	plt.figure(3)
	plt.contourf(x,y,EKE_HP,15,cmap=plt.cm.jet)
	plt.colorbar()
	
	plt.show()
	
	
	EKE_NEW = EKE_LP
	#mask1 = (kGrid*kGrid + lGrid*lGrid >= cutoffs[iFilter-1]*cutoffs[iFilter-1]) 
	#mask2 =  kGrid*kGrid + lGrid*lGrid < cutoffs[iFilter]*cutoffs[iFilter]
	
	#EKE_BP_FFT = fft_EKE.copy()
	#EKE_BP_FFT[:,mask1] = 0.0
	#EKE_BP_FFT[:,mask2] = 0.0
	
	#EKE_BP_FFT[mask1] = 0.0
	#EKE_BP_FFT[mask2] = 0.0
	
	
	#EKE_Int_FromFFT.append((EKE_BP_FFT*EKE_BP_FFT.conj()).sum())
	#
	#EKE_BP = ifftn(fftshift(EKE_BP_FFT),axes=[1,0])
	EKE_Integrated.append(np.sum(((EKE_HP[lowerIndex:upperIndex,bufferLength/2::]*EKE_HP[lowerIndex:upperIndex,bufferLength/2::].conj()))  ))
	
#	plt.figure(1)
	#plt.contourf(xWaveNumbers,yWaveNumbers,np.log10(fft_EKE*fft_EKE.conj()),15,cmap=plt.cm.jet)
#	plt.colorbar()
#	plt.figure(2)
#	plt.contourf(xWaveNumbers,yWaveNumbers,np.log10(EKE_BP_FFT*EKE_BP_FFT.conj()),15,cmap=plt.cm.jet)
#	plt.colorbar()
	
	#plt.figure(3)
	#plt.plot(xWaveNumbers[nX/2::],(fft_EKE[(nY/2)+(10),nX/2::])*fft_EKE[(nY/2)+(10),nX/2::].conj())
	#plt.figure(4)
	#plt.plot(yWaveNumbers[nY/2::],(fft_EKE[nY/2::,(nX/2)+10])*(fft_EKE[nY/2::,(nX/2)+10]).conj())
	
	
	#plt.figure(3)
	#plt.contourf(x,y,EHE_HP,15,cmap=plt.cm.jet)
	#plt.colorbar()
	
	plt.show()
		
	#EKE_NEW = EKE_LP
'''
#for iY in range(0,yWaveNumbers.shape[0]):
#	for iX in range(0,xWaveNumbers.shape[0]):
#		if not mask[iY,iX]:
#			print mask[iY,iX], iY, iX




print cutoffsX
for iFilter in range(0,len(cutoffs)):
	print iFilter
	#cutoffsX.append(cutOffWavelengthX)
	#cutoffsY.append(cutOffWavelengthY)
	cutOffkX  = cutoffs[iFilter]

	#cutOffkY = 1.0/cutOffWavelengthY
	
	print "Cut off wavelength", 1.0/cutOffkX
	print "Cut off wavenumber", cutOffkX
	
	#Compute kernal length of the filter

	#print nyquistFrequency
	tapWeightsX = firwin(kernalLength, cutOffkX, window='blackman', nyq=nyquistFrequencyX)
	tapWeightsY = firwin(kernalLength, cutOffkX, window='blackman', nyq=nyquistFrequencyY)
	
	EKE_LP_X = filtfilt(tapWeightsX, [1.0], EKE_NEW_X ,axis=1)
	EKE_LP_X = filtfilt(tapWeightsY, [1.0], EKE_LP_X ,axis=0)

	EHE_HP_X = EKE_NEW_X-EKE_LP_X
	#EHE_HP_Y = EKE_NEW_Y-EKE_LP_Y
	EKE_IntegratedX.append(np.sum(   (EHE_HP_X[lowerIndex::,bufferLength::]*EHE_HP_X[lowerIndex::,bufferLength::].conj())  ))
	#EKE_IntegratedY.append(np.sum(   (EHE_HP_Y[0,lowerIndex::,bufferLength::]*EHE_HP_Y[0,lowerIndex::,bufferLength::].conj())  ))

	print "Int EKE: ", EKE_IntegratedX[-1]
	
	#cutOffWavelengthX = 2.0 * cutOffWavelengthX 
	#cutOffWavelengthY = 2.0 * cutOffWavelengthY
	
	
	EKE_NEW_X = EKE_LP_X
	#EKE_NEW_Y = EKE_LP_Y

print cutoffsX
fig = plt.figure(1)
ax  = fig.add_subplot(1,1,1)
ax.contourf(x[kernalLength::],y[kernalLength::],EKE[:,kernalLength::,kernalLength::].mean(axis=0),15,cmap=plt.cm.jet)

#fig = plt.figure(2)
#ax  = fig.add_subplot(1,1,1)
#ax.contourf(x[kernalLength::],y[kernalLength::],EKE_LP_X[:,kernalLength::,kernalLength::].mean(axis=0),15,cmap=plt.cm.jet)

fig = plt.figure(3)
ax  = fig.add_subplot(1,1,1)
ax.plot(1.0/np.asarray(cutoffs), densityLayers[1]*np.asarray(EKE_IntegratedX)*float(nX*nY)*deltaX*deltaY*1.0e3,'k')

#fig = plt.figure(4)
#ax  = fig.add_subplot(1,1,1)
#ax.plot(1.0/np.asarray(cutoffs[10:-1]), EKE_Int_FromFFT,'k')

#ax.plot(cutoffsY, np.asarray(EKE_IntegratedY)*float(nX*nY)*deltaX*deltaY*1.0e3,'k--')
#ax.legend(('x-direction', ion'),
#           'upper right', shadow=True)
           
ax.set_ylabel('$\int\;$EKE (J.m$^{-2}$)')
ax.set_xlabel('Length Scale (km)')

plt.show()

dataFile.close()
