import numpy as np
import os
from scipy.io import netcdf as scipy_netcdf
from matplotlib import pyplot as plt
from matplotlib.path import Path
from scipy.signal import firwin, filtfilt, lfilter, find_peaks_cwt
from scipy.fftpack import fft, fftfreq, fftshift, ifft, fft2, ifft2
import PeakDetection as PD
from scipy.interpolate import griddata
from ContourMinimumSmoothing import SmoothContour
import matplotlib.gridspec as gridspec

def FindPeaksSpectrum(inputFFT):
	
	nPoints = inputFFT.shape[0]  
	smoothWindow = 5
	baselineWindow = 5
	smoothOrder = 2
	smoothFFT = PD.savitzky_golay(inputFFT, smoothWindow, smoothOrder)
	#peakind = find_peaks_cwt(inputFFT, np.arange(1,25))
	splitFFT = np.array_split(smoothFFT,int(nPoints/baselineWindow))
	splitAxis = np.array_split(np.arange(0,nPoints,1),int(nPoints/baselineWindow))
	
	
	baseLineValues  = np.zeros(len(splitFFT),dtype='float64')
	intepGridPoints = np.zeros(len(splitFFT),dtype='float64')

	for iSegment in range(0,len(splitFFT)):
		baseLineValues[iSegment] = np.mean(splitFFT[iSegment])
		intepGridPoints[iSegment] = splitAxis[iSegment][int(splitAxis[iSegment].shape[0]/2)]
	
	baseline = griddata(intepGridPoints, baseLineValues, np.arange(0,nPoints,1), method='linear')
	
	fftSmoothNoBaseline = smoothFFT - baseline
	firstDifference = np.diff(fftSmoothNoBaseline)
	maxima = np.diff((firstDifference > 0).view(np.int8))
	maxIndicies = np.concatenate((([0],) if firstDifference[0] < 0 else ()) +
                          (np.where(maxima == -1)[0] + 1,) +
                          (([len(fftSmoothNoBaseline)-1],) if firstDifference[-1] > 0 else ()))
	
	
	thresholdIndicies = []
	
	thresholdValue = 2.0*inputFFT.std()
	for iIndex in range(0,maxIndicies.shape[0]):
		if inputFFT[maxIndicies[iIndex]] > thresholdValue and maxIndicies[iIndex] != 0:
			print "thresholded: ", maxIndicies[iIndex]
			thresholdIndicies.append(maxIndicies[iIndex])
	if len(thresholdIndicies) > 1:
		thresholdIndicies.pop()
	
	
	#print len()
	#thresholdIndicies = np.asarray(thresholdIndicies).view(np.int8)
	#print "thresholding: ", thresholdIndicies
	'''
	plt.figure(1)
	plt.plot(inputFFT,'b')
	plt.plot(smoothFFT,'r')
	plt.plot(baseline,'g')
	
	plt.figure(2)
	plt.plot(smoothFFT - baseline)
	
	plt.figure(3)
	plt.plot((firstDifference))
	
	plt.figure(4)
	plt.plot(inputFFT,'b')
	plt.plot(2.0*inputFFT.std()*np.ones(nPoints,dtype='float64'),'r')
	plt.plot(thresholdIndicies,inputFFT[thresholdIndicies],'b*')
	
	#plt.plot(peakind,inputFFT[peakind],'b*')
	plt.show()
	'''
	#peakind = 1
	return thresholdIndicies
	'''
	nPoints = inputFFT.shape[-1]
	#inputSpectrum = inputFFT * inputFFT.conj()
	firstDeriv = np.diff(inputFFT*inputFFT)
	zero_crossings = np.where(np.diff(np.sign(firstDeriv)))[0]
	zero_crossing_thresh = []
	threshold = 3.0*np.std(np.abs(inputFFT))
	
	for iCrossing in range(0,len(zero_crossings)):
		if np.abs(inputFFT[zero_crossings[iCrossing]]) > threshold:
			zero_crossing_thresh.append(zero_crossings[iCrossing]+1)
	
	return zero_crossing_thresh
	'''
def create_empty_array_of_shape(shape):
    if shape: return [create_empty_array_of_shape(shape[1:]) for i in xrange(shape[0])]

baseDir = '/home/chris/GOLD/IsoTopo5/'
outputDirBaseString = 'output'
startYear = 40
endYear   = 60

topoFileName = '/home/chris/GOLD/IsoTopo_5_Processing/isoTopo_5.nc'
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

time        = np.zeros(nT_total,dtype='float64')
contourPositions = [ 0 ] * nT_total

iLayer = 4
fileCounter = 0

interpContourPositions = np.zeros([nT_total,nX_h],dtype='float64')
contourVal = 290.0


wavenumbers = fftshift(fftfreq(xh.shape[0],d=(xh[1]-xh[0])))


for iFileName in uvFilesList:

	nPointsToSmooth = 5

	fileObject = scipy_netcdf.netcdf_file(iFileName, 'r')
	time[fileCounter*nT_singleFile:(fileCounter+1)*nT_singleFile] = fileObject.variables['Time'][:]
	e = fileObject.variables['e'][:,:,:,:]
	u = fileObject.variables['u'][:,:,:,:]
	fileObject.close()

	eTimeMean = e.mean(axis=0)
	uTimeMean = u.mean(axis=0)
	
	fig = plt.figure(0)
	cs = plt.contour(xh,yh,-eTimeMean[iLayer,:,:],[contourVal])
	paths = cs.collections[0].get_paths()
	plt.close(fig)
	
	v = paths[0].vertices
	xPoints = v[:,0]
	yPoints = v[:,1] 
	
		
	xPointsRev = xPoints[::-1]
	yPointsRev = yPoints[::-1]
	meanContourPosition = np.vstack([xPointsRev,yPointsRev]).transpose()

	for iT in range(0,nT_singleFile):
	
		fig = plt.figure(0)
		cs = plt.contour(xh,yh,-e[iT,iLayer,:,:],[contourVal])
		paths = cs.collections[0].get_paths()
		plt.close(fig)
	
		maxPathLength = 0
		pathIndex     = 0
	
		for iPath in range(0,len(paths)):
			if (paths[iPath].vertices).shape[0] > maxPathLength:
				maxPathLength = (paths[iPath].vertices).shape[0]
				pathIndex = iPath
		v = paths[pathIndex].vertices
		xPoints = v[:,0]
		yPoints = v[:,1] 
	
		if xPoints[0] > xPoints[-1]: 
		
			xPointsRev = xPoints[::-1]
			yPointsRev = yPoints[::-1]
			contourPositions = np.vstack([xPointsRev,yPointsRev]).transpose()
		else:
			contourPositions = np.vstack([xPoints,yPoints]).transpose()
			
		
		interpContourPositions[iT + fileCounter*nT_singleFile,:] =  SmoothContour(contourPositions, xh)
		
	fileCounter = fileCounter + 1
'''
yLowerIndex = np.nonzero(yh>=250)[0][0]
yUpperIndex = np.nonzero(yh>1250)[0][0]


timeToPlot =nT_total-20 
valuesForPlot = np.arange(-1150,100,50)

meanDepth = 4000.0
topoFile = scipy_netcdf.netcdf_file(topoFileName, 'r')
topo        = meanDepth - topoFile.variables['topo'][:,:]

topoFile.close()

f, axarr = plt.subplots(2, sharex=True)
cs = axarr[0].contourf(xh,yh[yLowerIndex:yUpperIndex],e[-20,iLayer,yLowerIndex:yUpperIndex,:],valuesForPlot,cmap=plt.cm.jet)
axarr[0].contour(xh,yh[yLowerIndex:yUpperIndex],topo[yLowerIndex:yUpperIndex,:],10, colors = 'k', linestyle='dashed')
axarr[0].plot(xh,interpContourPositions[timeToPlot,:],linewidth=2.0,color='k')
axarr[0].yaxis.set_ticks([250,500,750,1000])
axarr[0].set_ylabel('y (km)')
axarr[0].annotate('a',xy=(0.01,1.025),xycoords='axes fraction',fontsize=15)

valuesForCbar = np.arange(-1000,0,250)
#fig.colorbar(cs,ax=axarr[0],shrink=.225, pad=.05, aspect=10,ticks=valuesForCbar)
cbar_ax = f.add_axes([0.91, 0.54, 0.015, 0.36])
cbar = f.colorbar(cs, cax=cbar_ax,ticks=valuesForCbar)
cbar.ax.set_title('SSH (m)')
#ax  = fig.add_subplot(2,1,2)
axarr[1].plot(xh,interpContourPositions[-24,:],'k')
axarr[1].plot(xh,interpContourPositions[-22,:]+500,'k')
axarr[1].plot(xh,interpContourPositions[-20,:]+1000,'k')
axarr[1].set_xlim(xh[0],xh[-1])
axarr[1].get_yaxis().set_visible(False)
axarr[1].set_xlabel('x (km)')
axarr[1].annotate('b',xy=(0.01,1.025),xycoords='axes fraction',fontsize=15)
timeString = 't=' + str(time[-24])
axarr[1].annotate(timeString,xy=(100,625),xycoords='data',fontsize=15)
timeString = 't=' + str(time[-22])
axarr[1].annotate(timeString,xy=(100,1125),xycoords='data',fontsize=15)
timeString = 't=' + str(time[-20])
axarr[1].annotate(timeString,xy=(100,1675),xycoords='data',fontsize=15)




#fig = plt.figure(2)
#ax  = fig.add_subplot(1,1,1,aspect = 'equal')
#ax.contourf(xh,yh[yLowerIndex:yUpperIndex],u[-20,iLayer,yLowerIndex:yUpperIndex,:],15,cmap=plt.cm.RdBu_r)
#ax.plot(xh,interpContourPositions[timeToPlot,:],linewidth=2.0,color='k')
plt.show()

'''
np.save("IsoTopo4_ContourPositions.npy", interpContourPositions)
#nterpContourPositions = np.load("IsoTopo4_ContourPositions.npy")
#print interpContourPositions.shape


#Perfom FFTs
print interpContourPositions[:,nX_h/2::].shape
contourFFT  = fftshift(fft2(interpContourPositions[:,nX_h/2::]))
print contourFFT.shape
#contourFFT = np.abs(contourFFT)
contourFFT  = (contourFFT*contourFFT.conj())
wavenumbers = fftshift(fftfreq(xh[nX_h/2::].shape[0],d=(xh[1]-xh[0])))
freqs       = fftshift(fftfreq(nT_total,d=(time[1]-time[0])))

nK = wavenumbers.shape[0]

#Compute the propagating spectrum
flipedFFT = contourFFT[:,::-1]

propSpectrum = contourFFT[nT_total/2::,nK/2:nK] - flipedFFT[nT_total/2::,nK/2:nK]

#Perform the frequency band analysis

nFreqBands = 100

 
splitFreqArray = np.array_split(freqs[nT_total/2::],nFreqBands+1)

bandLimitedSpectra1 =  np.zeros([nFreqBands+1,nK/2],dtype='float64')
bandLimitedSpectra2 =  np.zeros([nFreqBands+1,nK/2],dtype='float64')

spectralPeakTotal = []
spectralPeakProp  = []
wavenumbersRight = wavenumbers[nK/2:nK]

print freqs[nT_total/2::]
for iBand in range(0,nFreqBands+1):

	lowerIndex = np.nonzero(freqs[nT_total/2::]>=splitFreqArray[iBand][0])[0][0]
	upperIndex = np.nonzero(freqs[nT_total/2::]>=splitFreqArray[iBand][-1])[0][0]
	
	
	print splitFreqArray[iBand][0]
	print splitFreqArray[iBand][-1]
	
	print nT_total/2
	print lowerIndex
	print upperIndex
	print (contourFFT[nT_total/2+lowerIndex:nT_total/2+upperIndex,:].mean(axis=0)).shape
	
	bandLimitedSpectra1[iBand,:] = contourFFT[nT_total/2+lowerIndex:nT_total/2+upperIndex,nK/2:nK].mean(axis=0)
	bandLimitedSpectra2[iBand,:] = (propSpectrum[lowerIndex:upperIndex,:]).mean(axis=0)
	
	#spectralPeakTotal.append(np.argmax(bandLimitedSpectra1[iBand,:]))
	spectralPeakTotal.append(FindPeaksSpectrum(bandLimitedSpectra1[iBand,:]))
	spectralPeakProp.append(FindPeaksSpectrum(bandLimitedSpectra2[iBand,:]))
#	
#	plt.figure(1)
#	plt.plot(wavenumbersRight, bandLimitedSpectra1[iBand,:],'b')
#	plt.plot(wavenumbersRight[spectralPeakTotal[-1]],bandLimitedSpectra1[iBand,spectralPeakTotal[-1]],'b*')
#	plt.show()
	
	
fig = plt.figure(0)
#ax  = fig.add_subplot(1,1,1)
gs = gridspec.GridSpec(2, 1,height_ratios=[2,1])
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])
#ax3 = ax1.twiny()



cs = ax1.contourf(wavenumbersRight[0:wavenumbersRight.shape[0]/2],freqs[nT_total/2::],np.log10(contourFFT[nT_total/2::,nK/2:3*nK/4]),15,cmap=plt.cm.jet)
ax1.plot(wavenumbersRight[0:wavenumbersRight.shape[0]/2],splitFreqArray[2][0]*np.ones( wavenumbersRight.shape[0]/2,dtype='float64'),'k')
ax1.plot(wavenumbersRight[0:wavenumbersRight.shape[0]/2],splitFreqArray[2][-1]*np.ones(wavenumbersRight.shape[0]/2,dtype='float64'),'k')
ax1.set_xlim(wavenumbersRight[0],wavenumbersRight[(wavenumbersRight.shape[0]/2)-1])
ax1.set_ylim(freqs[nT_total/2],freqs[-1])

xTicksWavelength = np.asarray([1000.0,500.0,250.0,100.0,50.0,25.0])

#ax1.set_xticklabels(xTicksWavelength)
#ax2.set_xticks(1.0/xTicksWavelength)
ax1.set_ylabel("Period  (days)")
ax1.annotate('a',xy=(-0.05,1.025),xycoords='axes fraction',fontsize=15)
#ax3.set_xticks(1.0/xTicksWavelength)
#ax3.set_xticklabels(xTicksWavelength)
#ax3.set_xlabel(r"Wavelength (km)")
#ax1.set_xticks(np.arange(0.01,0.))
ax1.get_xaxis().set_visible(False)
ax1.set_yticks(np.arange(0.01,0.041,0.01))
ax1.set_yticklabels((1.0/np.arange(0.01,0.041,0.01)).astype(np.int16))


#fig.colorbar(cs,ax=ax1)
ax2.plot(wavenumbersRight[0:wavenumbersRight.shape[0]/2],bandLimitedSpectra1[2][0:bandLimitedSpectra1[2].shape[0]/2],'k')
ax2.set_xlim(wavenumbersRight[0],wavenumbersRight[(wavenumbersRight.shape[0]/2)-1])
ax2.set_xlabel('Wavelength (km)')
ax2.set_ylabel('Ave Spectral Density')
ax2.annotate('b',xy=(-0.05,1.1),xycoords='axes fraction',fontsize=15)
ax2.set_xticks(np.hstack((np.asarray(1.0e-3),np.arange(0.005,0.036,0.005))))
ax2.set_xticklabels((1.0/np.hstack((np.asarray(1.0e-3),np.arange(0.005,0.036,0.005)))).astype(np.int16))


fig = plt.figure(1)
ax  = fig.add_subplot(1,1,1)
cs = ax.contourf(wavenumbersRight,freqs[nT_total/2::],np.log10(propSpectrum),15,cmap=plt.cm.jet)
fig.colorbar(cs,ax=ax)

bandToPlot = 2
print "lower band freq: ", splitFreqArray[iBand][0]
print "upper band freq: ", splitFreqArray[iBand][-1]

print "lower band period: ", 1.0/splitFreqArray[iBand][0]
print "upper band period: ", 1.0/splitFreqArray[iBand][-1]

'''
fig = plt.figure(2)
for iBand in range(1,nFreqBands+1):
	ax  = fig.add_subplot(nFreqBands,1,iBand)
	print iBand
	print iBand-1
	ax.plot(wavenumbers[nK/2:nK],(bandLimitedSpectra1[iBand,:]),'b')
	ax.plot(wavenumbersRight[spectralPeakTotal[iBand]],bandLimitedSpectra1[iBand,spectralPeakTotal[iBand]],'b*')

'''
print spectralPeakTotal[0]

fig = plt.figure(3)
ax = fig.add_subplot(1,1,1)
for iBand in range(1,nFreqBands+1):
	revBand = nFreqBands - iBand + 1
	aveFreq = 0.5*(splitFreqArray[iBand][0]+splitFreqArray[iBand][-1])
	wavesToGet = wavenumbersRight[np.asarray(spectralPeakTotal[iBand]).astype(np.int8)]
	ax.plot(-1.0/wavenumbersRight[np.asarray(spectralPeakTotal[iBand]).astype(np.int8)], aveFreq/(wavesToGet),'ko')
	
	ax.set_xticks(np.arange(-900,1,100))
	ax.set_xticklabels(np.arange(0,901,100))
	ax.set_xlabel("wavelength (km)")
	ax.set_ylabel("Phase speed (km/day)")
	
	#ax.plot(wavenumbersRight[np.asarray(spectralPeakTotal[iBand])], aveFreq,'ko')

fig = plt.figure(4)
ax = fig.add_subplot(1,1,1)
for iBand in range(1,nFreqBands+1):
	revBand = nFreqBands - iBand + 1
	aveFreq = 0.5*(splitFreqArray[iBand][0]+splitFreqArray[iBand][-1])
	wavesToGet = wavenumbersRight[np.asarray(spectralPeakProp[iBand]).astype(np.int8)]
	
	ax.plot(-1.0/wavenumbersRight[np.asarray(spectralPeakProp[iBand]).astype(np.int8)], aveFreq/wavesToGet,'ko')
	#ax.plot(wavenumbersRight[np.asarray(spectralPeakTotal[iBand])], aveFreq,'ko')
	
	ax.set_xticks(np.arange(-900,1,100))
	ax.set_xticklabels(np.arange(0,901,100))
	ax.set_xlabel("wavelength (km)")
	ax.set_ylabel("Phase speed (km/day)")
plt.show()






