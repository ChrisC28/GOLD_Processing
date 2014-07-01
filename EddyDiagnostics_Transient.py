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


baseDir = '/home/chris/GOLD/IsoTopo4/'
outputDirBaseString = 'output'
startYear = 40
endYear   = 42


outputDirs = []

for iYear in range(startYear, endYear+1):
	outputDirs.append(outputDirBaseString + "%03d" % iYear + '/')
print outputDirs
fileStartingString = 'ave_prog__'
progFileStartingString = 'prog__'

uvFilesList = []
pvFilesList = []
for iYear in range(len(outputDirs)):
	#print baseDir + outputDirs[iYear]
	for file in os.listdir(baseDir + outputDirs[iYear]):
		if os.path.isfile(baseDir + outputDirs[iYear] + file) and file.startswith(fileStartingString):
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

print nT_singleFile
print nT_total
#print uvFilesList
tileSizeX = 10
tileSizeY = 10

n_Xslice   = 32
xqSliceSize = nX_q / n_Xslice
xhSliceSize = nX_h / n_Xslice
 


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





time = np.zeros(nT_total,dtype='float64')
u    = np.zeros((nT_total,nZ_l,3,xqSliceSize+2),dtype='float64')
v    = np.zeros((nT_total,nZ_l,3,xhSliceSize+2),dtype='float64')
e    = np.zeros((nT_total,3,xhSliceSize+2),dtype='float64')
h    = np.zeros((nT_total,nZ_l,3,xhSliceSize+2),dtype='float64')
q    = np.zeros((nT_total,nZ_l,3,xhSliceSize+2),dtype='float64')


uLowPass  = np.zeros((nT_total,nZ_l/2,xqSliceSize),dtype='float64')
vLowPass  = np.zeros((nT_total,nZ_l/2,xhSliceSize),dtype='float64')

uTransientEddy = np.zeros((nT_total,nZ_l/2,xqSliceSize),dtype='float64')
vTransientEddy = np.zeros((nT_total,nZ_l/2,xhSliceSize),dtype='float64')
eTransientEddy = np.zeros((nT_total,nZ_i/2+1,xqSliceSize),dtype='float64')


uTimeMean = np.zeros((nZ_l,nY_h-2,nX_q),dtype='float64')
vTimeMean = np.zeros((nZ_l,nY_q-2,nX_h),dtype='float64')
hTimeMean = np.zeros((nZ_l,nY_q-2,nX_h),dtype='float64')
mTimeMean = np.zeros((nZ_l,nY_q-2,nX_h),dtype='float64')
qTimeMean = np.zeros((nZ_l,nY_q-2,nX_h),dtype='float64')


reynoldsStressTimeMean = np.zeros((nZ_l,nY_q-2,nX_h),dtype='float64')
potFluxTimeMean        = np.zeros((nZ_l,nY_q-2,nX_h),dtype='float64')


for iY in range(1,nY_h-1):
	print float(iY)/float(nY_h-1) 	
	for iX in range(0,n_Xslice):

		fileCounter  = 0
		for iFile in range(0,len(uvFilesList)):
			
			fileObject = scipy_netcdf.netcdf_file(uvFilesList[iFile], 'r')
			startXIndex = iX*xqSliceSize
			if iX ==0:
				#endXIndex = ((iX+1) * xqSliceSize)+1
				xSliceIndicies = range(startXIndex,((iX+1) * xqSliceSize)+1,1)
				xSliceIndicies.insert(0,nX_q-1)
				
			elif (iX != 0) and (iX < n_Xslice-1):
				#endXIndex = ((iX+1) * xqSliceSize)+1
				xSliceIndicies = range(startXIndex-1,((iX+1) * xqSliceSize)+1,1)
				
			else:
				endXIndex = ((iX+1) * xqSliceSize)
				xSliceIndicies = range(startXIndex-1,((iX+1) * xqSliceSize),1)
				xSliceIndicies.append(0)
			endXIndex = (iX+1)*xqSliceSize
		
			u[fileCounter*nT_singleFile:(fileCounter+1)*nT_singleFile,:,:,:]  = fileObject.variables['u'][:,:,iY-1:iY+2,xSliceIndicies]
			v[fileCounter*nT_singleFile:(fileCounter+1)*nT_singleFile,:,:,:]  = fileObject.variables['v'][:,:,iY-1:iY+2,xSliceIndicies]
			eTemp = fileObject.variables['e'][:,:,iY-1:iY+2,xSliceIndicies]
			e[fileCounter*nT_singleFile:(fileCounter+1)*nT_singleFile,:,:]  = eTemp[:,0,:,:]
			h[fileCounter*nT_singleFile:(fileCounter+1)*nT_singleFile,:,:,:]  = fileObject.variables['h'][:,:,iY-1:iY+2,xSliceIndicies]
			
			fileObject.close()
			
			fileObject = scipy_netcdf.netcdf_file(pvFilesList[iFile], 'r')
			q[fileCounter*nT_singleFile:(fileCounter+1)*nT_singleFile,:,:,:]  = fileObject.variables['PV'][:,:,iY-1:iY+2,xSliceIndicies]
			
			fileObject.close()
			
			fileCounter = fileCounter + 1
			
		
		for iLayer in range(0,nZ_l):
			h[:,iLayer,:,:] = h[:,iLayer,:,:] - refThickness[iLayer]
			
		#Interpolate to 'h' grid nodes
		#u_h = 0.5*(u[:,:,1::]+u[:,:,0:-1])
		#v_h = 0.5*(np.squeeze(v[:,:,1::,:]+v[:,:,0:-1,:]))
	
		#print "filtering"
		uLowPass = lfilter(tapWeights,[1.0],u, axis=0)	
		vLowPass = lfilter(tapWeights,[1.0],v, axis=0)	
		eLowPass = lfilter(tapWeights,[1.0],e, axis=0)	
		hLowPass = lfilter(tapWeights,[1.0],h, axis=0)	
		qLowPass = lfilter(tapWeights,[1.0],q, axis=0)	
		#print "finished filtering"
	
		uTransientEddy = u - uLowPass
		vTransientEddy = v - vLowPass
		eTransientEddy = e - eLowPass
		hTransientEddy = h - hLowPass
		qTransientEddy = q - qLowPass
		
		montPotEddy = MP.ComputeMontPot(eTransientEddy,hTransientEddy,reducedGravity,densityLayers, RHO0,G0)
		
		#print vTimeMean[:,iY-1,iX*xqSliceSize:(iX+1)*xqSliceSize].shape
		#print np.squeeze(v_h.mean(axis=0)).shape
		uTimeMean[:,iY-1,iX*xqSliceSize:(iX+1)*xqSliceSize] = np.squeeze((u[:,:,1,1:u.shape[-1]-1]).mean(axis=0))
		vTimeMean[:,iY-1,iX*xqSliceSize:(iX+1)*xqSliceSize] = np.squeeze((v[:,:,1,1:u.shape[-1]-1]).mean(axis=0))
		hTimeMean[:,iY-1,iX*xqSliceSize:(iX+1)*xqSliceSize] = np.squeeze((h[:,:,1,1:u.shape[-1]-1]).mean(axis=0))
		mTimeMean[:,iY-1,iX*xqSliceSize:(iX+1)*xqSliceSize] = np.squeeze((montPotEddy[:,:,1,1:u.shape[-1]-1]).mean(axis=0))
		qTimeMean[:,iY-1,iX*xqSliceSize:(iX+1)*xqSliceSize] = np.squeeze((qTransientEddy[:,:,1,1:q.shape[-1]-1]).mean(axis=0))
		
		
		#EP_Vector = MP.ComputeEPvector(uTransientEddy,vTransientEddy,montPotEddy, uLowPass[:,:,1,:], vLowPass[:,:,1,:], 
		#				   refThickness, densityLayers, deltaX*1.0e3, deltaY*1.0e3, 10.0*G0,F0,RHO0)
		
		#EP_X_TimeMean[:,iY-1,iX*xqSliceSize:(iX+1)*xqSliceSize] = EP_Vector[0].mean(axis=0)
		#EP_Y_TimeMean[:,iY-1,iX*xqSliceSize:(iX+1)*xqSliceSize] = EP_Vector[1].mean(axis=0)
		#EP_Z_TimeMean[:,iY-1,iX*xqSliceSize:(iX+1)*xqSliceSize] = EP_Vector[2].mean(axis=0)
		
		
		#EKE_Transient_TimeMean[:,iY-1,iX*xqSliceSize:(iX+1)*xqSliceSize] =  0.5*( uTransientEddy[:,:,1,1:u.shape[-1]-1]*uTransientEddy[:,:,1,1:u.shape[-1]-1] + 
														#						  vTransientEddy[:,:,1,1:v.shape[-1]-1]*vTransientEddy[:,:,1,1:v.shape[-1]-1] ).mean(axis=0)
		#PE_Transient_TimeMean[:,iY-1,iX*xqSliceSize:(iX+1)*xqSliceSize]  = 0.5*( hTransientEddy[:,:,1,1:h.shape[-1]-1]*hTransientEddy[:,:,1,1:h.shape[-1]-1] ).mean(axis=0)
		reynoldsStress, potFlux = MP.ComputeEP_Components(uTransientEddy,vTransientEddy,montPotEddy, uLowPass[:,:,1,:], vLowPass[:,:,1,:], 
						   refThickness, densityLayers, deltaX*1.0e3, deltaY*1.0e3, 10.0*G0,F0,RHO0)

		reynoldsStressTimeMean[:,iY-1,iX*xqSliceSize:(iX+1)*xqSliceSize] = reynoldsStress.mean(axis=0)
		potFluxTimeMean[:,iY-1,iX*xqSliceSize:(iX+1)*xqSliceSize]        = potFlux.mean(axis=0)

#np.save('EP_X_Transient_IsoTopo4.npy',EP_X_TimeMean)
#np.save('EP_Y_Transient_IsoTopo4.npy',EP_Y_TimeMean)
#np.save('EP_Z_Transient_IsoTopo4.npy',EP_Z_TimeMean)

#np.save('EKE_Transient_IsoTopo4.npy',EKE_Transient_TimeMean)
#np.save('PE_Transient_IsoTopo4.npy',PE_Transient_TimeMean)
#np.save('Q_Transient_IsoTopo4.npy',qTimeMean)



postProcessedOutput = scipy_netcdf.netcdf_file('waveActivity_Components_Iso4.nc', 'w')

xDim = postProcessedOutput.createDimension('x', xh.shape[0])
yDim = postProcessedOutput.createDimension('y', yh.shape[0]-2)
zlDim = postProcessedOutput.createDimension('zl', zl.shape[0])



xVar = postProcessedOutput.createVariable('x','f4',('x',))
yVar = postProcessedOutput.createVariable('y','f4',('y',))
zlVar = postProcessedOutput.createVariable('zl','f4',('zl',))

xVar[:] = xh
yVar[:] = yh[1:nY_h-1]
zlVar[:] = zl

RS_X_Var = postProcessedOutput.createVariable('ReynoldsStress','f4',('zl','y','x'))
PFlux_X_Var = postProcessedOutput.createVariable('PotentialFlux','f4',('zl','y','x'))

#EPXVar = postProcessedOutput.createVariable('EP_X','f4',('zl','y','x'))
#EPYVar = postProcessedOutput.createVariable('EP_Y','f4',('zl','y','x'))
#EPZVar = postProcessedOutput.createVariable('EP_Z','f4',('zl','y','x'))

#PE_Var = postProcessedOutput.createVariable('PE','f4',('zl','y','x'))



RS_X_Var[:,:,:] = reynoldsStressTimeMean
PFlux_X_Var[:,:,:] = potFluxTimeMean


postProcessedOutput.close()

	
	
layerToPlot = 1
fig = plt.figure(1)
ax  = fig.add_subplot(1,1,1)
cs = ax.contourf(xq,yh[1:nY_h-1],uTimeMean[layerToPlot,:,:],15,cmap=plt.cm.jet)
fig.colorbar(cs,ax=ax)

              
fig = plt.figure(2)
ax  = fig.add_subplot(1,1,1)
cs = ax.contourf(xq,yh[1:nY_h-1],vTimeMean[layerToPlot,:,:],15,cmap=plt.cm.jet)
fig.colorbar(cs,ax=ax)

fig = plt.figure(3)
ax  = fig.add_subplot(1,1,1)
cs = ax.contourf(xq,yh[1:nY_h-1],EP_Z_TimeMean[layerToPlot,:,:],15,cmap=plt.cm.jet)
fig.colorbar(cs,ax=ax)


fig = plt.figure(4)
ax  = fig.add_subplot(1,1,1)
cs = ax.contourf(xq,yh[1:nY_h-1],EKE_Transient_TimeMean[layerToPlot,:,:],15,cmap=plt.cm.jet)
fig.colorbar(cs,ax=ax)

fig = plt.figure(5)
ax  = fig.add_subplot(1,1,1)
cs = ax.contourf(xq,yh[1:nY_h-1],PE_Transient_TimeMean[layerToPlot,:,:],15,cmap=plt.cm.jet)
fig.colorbar(cs,ax=ax)


fig = plt.figure(6)
ax  = fig.add_subplot(1,1,1)
cs = ax.contourf(xq,yq[1:nY_h-1],qTimeMean[layerToPlot,:,:],15,cmap=plt.cm.jet)
fig.colorbar(cs,ax=ax)


normEP = np.sqrt(EP_X_TimeMean**2+EP_Y_TimeMean**2)

maskArray = normEP[layerToPlot,:,:] < 1.0*normEP[layerToPlot,:,:].std()

EP_X_mask = np.ma.masked_array(EP_X_TimeMean[layerToPlot,:,:], mask=maskArray)
EP_Y_mask = np.ma.masked_array(EP_Y_TimeMean[layerToPlot,:,:], mask=maskArray)



fig = plt.figure(7)
ax = fig.add_subplot(1,1,1)
cs = ax.contourf(xh,yh[1:nY_h-1],EP_Z_TimeMean[layerToPlot+2,:,:],15,cmap=plt.cm.jet)
fig.colorbar(cs,ax=ax)

X,Y = np.meshgrid(xh,yh[1:nY_h-1])

#fig = plt.figure(4)
Q = ax.quiver( X[::20, ::20], Y[::20, ::20], EP_X_mask[::20, ::20], 
                                          EP_X_mask[::20, ::20], 
										  pivot='mid', color='k',units='width',scale=2.5)


plt.show()
