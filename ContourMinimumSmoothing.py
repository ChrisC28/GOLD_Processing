import numpy as np
import os
from scipy.io import netcdf as scipy_netcdf
from matplotlib import pyplot as plt
from matplotlib.path import Path
from scipy.signal import firwin, filtfilt, lfilter
from scipy.fftpack import fft, fftfreq, fftshift, ifft, fft2, ifft2
import PeakDetection as PD
from scipy.interpolate import griddata

def SmoothContour(contourPositions, xh):
	
	lastX = contourPositions[0,0]
	nextX = contourPositions[0,0]

	enterBendFlag = False
	midBendFlag  = False

	startBendIndicies = []
	endBendIndicies   = []
	midBendIndicies   = []
	startBendValue = -99.0
	
	minXValue = -10.0
	for iStep in range(1,contourPositions[:,0].shape[0]):
		nextX = contourPositions[iStep,0]
		
		if (nextX < lastX) and (not enterBendFlag):
		
			enterBendFlag = True
			midBendFlag   = True
			startBendIndicies.append(iStep)
			startBendValue = nextX
			minXValue   = nextX
			minXIndex = iStep 
		
		if (nextX < minXValue) and (enterBendFlag): 
			minXValue = nextX
			minXIndex = iStep 
		
		if (nextX > startBendValue) and (enterBendFlag):
		
			enterBendFlag = False
			midBendIndicies.append(minXIndex)
			endBendIndicies.append(iStep)
	
		lastX = nextX
	
	firstMinStartIndicies = []
	
	for iIndex in midBendIndicies:
	
		minXValue = contourPositions[iIndex,0]
		firstMinStartIndicies.append(np.nonzero(contourPositions[:,0]>=minXValue)[0][0] )
	
	newStartIndex = []
	
	for iSegment in range(0,len(startBendIndicies)):
		#print "iSegment: ", iSegment
		#print firstMinStartIndicies
		if not firstMinStartIndicies:
			newStartIndex.append(startBendIndicies[iSegment])
			#print "blah!"
			#plt.figure(1)
			#plt.plot(contourPositions[:,0],contourPositions[:,1])
			#plt.plot(contourPositions[startBendIndicies,0],contourPositions[startBendIndicies,1],'r*')
			#plt.show()
		
		else:
			minXValue = contourPositions[firstMinStartIndicies[iSegment],0]

			startXValue = contourPositions[startBendIndicies[iSegment],0]

			newStartXValue = 0.5*(minXValue+startXValue)
	
			newStartIndex.append(np.nonzero(contourPositions[firstMinStartIndicies[iSegment]:startBendIndicies[iSegment],0]>=newStartXValue)[0][0] 
								+ firstMinStartIndicies[iSegment] )

	newMidIndex = []

	for iSegment in range(0,len(startBendIndicies)):
	
		minXValue = contourPositions[newStartIndex[iSegment],0]
		minYValue = contourPositions[newStartIndex[iSegment],1]

	
		startXValue = contourPositions[endBendIndicies[iSegment],0]
		startYValue = contourPositions[endBendIndicies[iSegment],1]
	
		newStartXValue = 0.5*(minXValue+startXValue)
		newStartYValue = 0.5*(minYValue+startYValue)
		newMidIndex.append(np.nonzero(contourPositions[startBendIndicies[iSegment]:endBendIndicies[iSegment],0]<newStartXValue)[0][-1]  
	                       + startBendIndicies[iSegment] )

	
	startingIndex = 0
	for iSegment in range(0,len(startBendIndicies)):
	
		if 0 == iSegment: 
	
			newContour = np.vstack([contourPositions[startingIndex:newStartIndex[iSegment],0], 
			     					contourPositions[startingIndex:newStartIndex[iSegment],1]])
		else:
		
			contourSlice = np.vstack([contourPositions[startingIndex:newStartIndex[iSegment],0], 
			     					  contourPositions[startingIndex:newStartIndex[iSegment],1]])
			newContour = np.hstack([newContour,contourSlice])
	
		midPoint = np.vstack([contourPositions[newMidIndex[iSegment],0], 
	                          contourPositions[newMidIndex[iSegment],1]])
	
		newContour = np.hstack([newContour,midPoint])
	                        
		startingIndex = endBendIndicies[iSegment]                
	
	
	if len(startBendIndicies) != 0:
		contourSlice = np.vstack([contourPositions[endBendIndicies[-1]::,0], 
								  contourPositions[endBendIndicies[-1]::,1]])

		newContour = np.hstack([newContour,contourSlice])
		newContour = newContour.transpose()
	else:
		newContour = contourPositions
	
	
	return griddata(newContour[:,0], newContour[:,1], xh, method='linear')


