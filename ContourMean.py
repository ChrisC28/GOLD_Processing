import numpy as np
import os
from scipy.io import netcdf as scipy_netcdf
from matplotlib import pyplot as plt
from scipy.signal import lfilter, firwin, filtfilt, freqz, convolve
from scipy.interpolate import griddata

def create_empty_array_of_shape(shape):
    if shape: return [create_empty_array_of_shape(shape[1:]) for i in xrange(shape[0])]


NORTH_LAT = -35.0
SOUTH_LAT = -65

WEST_LONG = 0.0
EAST_LONG = 180.0

EARTH_RADIUS = 6380.0e3
DEG_2_RAD = np.pi/180.0
EARTH_ROTATION_RATE = 7.2921150e-5
G  = 9.81
oneOnEarthRad = 1.0/(EARTH_RADIUS)


startYear = 39
endYear   = 69


baseDir = "/home/chris/GOLD/IsoTopo5/"
outputDirBaseString = "output"
variableName = 'e'
uvFilesList = []
outputDirs = []

for iYear in range(startYear, endYear+1):
	outputDirs.append(outputDirBaseString + "%03d" % iYear + '/')
print outputDirs
fileStartingString = 'ave_prog__'

for iYear in range(len(outputDirs)):
	#print baseDir + outputDirs[iYear]
	for file in os.listdir(baseDir + outputDirs[iYear]):
		if os.path.isfile(baseDir + outputDirs[iYear] + file) and file.startswith(fileStartingString):
			uvFilesList.append(baseDir + outputDirs[iYear] + file)  
uvFilesList = sorted(uvFilesList)
print uvFilesList

numStepsToAverage = len(uvFilesList)


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

deltaX = xh[1]-xh[0]
deltaY = yh[1]-yh[0]


#sshFile = scipy_netcdf.netcdf_file(uvFilesList[0], 'r')
#lon = sshFile.variables['Longitude'][:]
#lat = sshFile.variables['Latitude'][:]
#iUpperLat = np.nonzero(lat > NORTH_LAT)[0][0]
#iLowerLat = np.nonzero(lat >= SOUTH_LAT)[0][0]

#iLeftLon  = np.nonzero(lon > WEST_LONG)[0][0]
#iRightLon = np.nonzero(lon >= EAST_LONG)[0][0]


#lon = lon[iLeftLon-1:iRightLon+1]
#lat = lat[iLowerLat:iUpperLat+1]

#nLat = lat.shape[0]
#nLon = lon.shape[0]
#deltaLat = lat[1]-lat[0]
#deltaLon = lon[1]-lon[0]
#print deltaLat
#print deltaLon

layer = 4
averagedSSH  = np.zeros([nY_h,nX_h],dtype='float64')
averagedSpeed  = np.zeros([nY_h,nX_h],dtype='float64')

for fileToGet in range(0,numStepsToAverage):
	print  uvFilesList[fileToGet]
	sshFile = scipy_netcdf.netcdf_file(uvFilesList[fileToGet], 'r')
	ssh = sshFile.variables[variableName][:,layer,:,:]
	averagedSSH = averagedSSH +  ssh.mean(axis=0)
	u = sshFile.variables['u'][:,layer,:,:]
	v = sshFile.variables['u'][:,layer,:,:]

	averagedSpeed = averagedSpeed +  (np.sqrt(u*u)).mean(axis=0)

	sshFile.close()
averagedSSH = averagedSSH/float(numStepsToAverage)
averagedSSH = np.ma.array(averagedSSH,mask=np.isnan(averagedSSH))
averagedSpeed = averagedSpeed/float(numStepsToAverage)
averagedSpeed = np.ma.array(averagedSpeed,mask=np.isnan(averagedSpeed))



'''
averagedOutputFileName = 'aveSSH_OFES.nc'

aveOutFile = scipy_netcdf.netcdf_file(averagedOutputFileName, 'w')
xDim = aveOutFile.createDimension('lon',length=nLon)
yDim = aveOutFile.createDimension('lat',length=nLat)

xVar = aveOutFile.createVariable('lon','f4',('lon',))
yVar = aveOutFile.createVariable('lat','f4',('lat',))
aveSSHVar = aveOutFile.createVariable('ssh','f4',('lat','lon'))

xVar[:] = lon
yVar[:] = lat
aveSSHVar[:,:] = averagedSSH


aveOutFile.close()
'''



firstDerivFDcoeffs  = np.asarray( [1.0, 0.0, -1.0] ) / (2.0)

#firstDerivFDcoeffs  = np.asarray( [1.0,-9.0,45.0, 0.0, -45.0,9.0,-1.0] ) / 60.0
d_SSH_dX = np.zeros([nY_h,nX_h],dtype='float64')
d_SSH_dY = np.zeros([nY_h,nX_h],dtype='float64')


for iY in range(0,nY_h):
	#oneOnCos = 1.0/(np.cos(DEG_2_RAD * lat[iY]))
	#D_SSH_Dlon[iY,:]  = oneOnEarthRad * oneOnCos * convolve(averagedSSH[iY,:], firstDerivFDcoeffs/(DEG_2_RAD*deltaLon), mode='same') * 100.0e3
	d_SSH_dX[iY,:]   = convolve(averagedSSH[iY,:], firstDerivFDcoeffs/(deltaX), mode='same') * 100.0e3


for iX in range(0,nX_h):

	d_SSH_dY[:,iX]   = oneOnEarthRad * convolve(averagedSSH[:,iX], firstDerivFDcoeffs/(deltaY), mode='same') * 100.0e3


gradMag = np.sqrt(d_SSH_dX*d_SSH_dX + d_SSH_dY*d_SSH_dY)




sshValues = np.arange(100.0,1000.0,10)
aveAlongSSHGradientValues = np.zeros(sshValues.shape[0],dtype='float64')
aveAlongSSHSpeedValues = np.zeros(sshValues.shape[0],dtype='float64')



fig=plt.figure(0)
cs = plt.contour(xh,yh,-averagedSSH,sshValues)
plt.close(fig)
contourPoints = []
grid_x, grid_y = np.meshgrid(xh, yh)
for iContour in range(0,len(cs.collections)):
	paths = cs.collections[iContour].get_paths()
	
	maxPathLength = 0
	pathIndex = 0
	
	for iPath in range(0,len(paths)):
		if (paths[iPath].vertices).shape[0] > maxPathLength:
			maxPathLength = (paths[iPath].vertices).shape[0]
			pathIndex = iPath
	#print len(paths)
	vPath = paths[pathIndex].vertices
	contourPoints.append(vPath)
#plt.figure(0)
	gradMag = np.ma.masked_array(gradMag,mask=np.isnan(gradMag))
	
	points = np.vstack([grid_x.reshape(nY_h*nX_h),grid_y.reshape(nY_h*nX_h)]).transpose()
	interpPoints = griddata(points, gradMag.reshape(nY_h*nX_h), vPath, method='nearest')
	
	aveAlongSSHGradientValues[iContour] =  interpPoints.mean()
	interpPoints = griddata(points, averagedSpeed.reshape(nY_h*nX_h), vPath, method='nearest')
	
	aveAlongSSHSpeedValues[iContour] =  interpPoints.mean()

vals = [290.0]
print sshValues.shape
print aveAlongSSHGradientValues.shape
plt.figure(0)
plt.plot(sshValues,aveAlongSSHSpeedValues)

plt.figure(1)
plt.plot(sshValues,aveAlongSSHGradientValues)


fig = plt.figure(2)
ax  = fig.add_subplot(1,1,1)
CS = ax.contourf(xh,yh,-ssh[0,:,:],15,cmap=plt.cm.jet)
ax.contour(xh,yh,-ssh[0,:,:],vals,colors='k')
cbar = plt.colorbar(CS)

fig = plt.figure(3)
ax  = fig.add_subplot(1,1,1)
CS = ax.contourf(xh,yh,u[0,:,:],15,cmap=plt.cm.jet)
ax.contour(xh,yh,-ssh[0,:,:],vals,colors='k')

cbar = plt.colorbar(CS)
plt.show()



exit()


frontSSHValues = [-1.35,-0.85,-0.27]

numTimeSteps = 20


frontPositions = create_empty_array_of_shape([numTimeSteps-1,len(frontSSHValues),1]) 


for fileToGet in range(1,numTimeSteps):
	#print  dataFilesList[fileToGet]
	sshFile1 = scipy_netcdf.netcdf_file(dataFilesList[(5*fileToGet)-1], 'r')
	sshFile2 = scipy_netcdf.netcdf_file(dataFilesList[(5*fileToGet)], 'r')
	sshFile3 = scipy_netcdf.netcdf_file(dataFilesList[(5*fileToGet)+1], 'r')

	timeAveSSH = (1.0/3.0)*(sshFile1.variables['eta'][0,iLowerLat:iUpperLat+1,iLeftLon-1:iRightLon+1] +
					        sshFile2.variables['eta'][0,iLowerLat:iUpperLat+1,iLeftLon-1:iRightLon+1] +
					        sshFile3.variables['eta'][0,iLowerLat:iUpperLat+1,iLeftLon-1:iRightLon+1]) 
	
	sshFile1.close()
	sshFile2.close()
	sshFile3.close()
	fig=plt.figure(0)
	cs = plt.contour(lon,lat,timeAveSSH,frontSSHValues)
	plt.close(fig)
	
	for iContour in range(0,len(cs.collections)):
		paths = cs.collections[iContour].get_paths()
	
		maxPathLength = 0
		pathIndex = 0
	
		for iPath in range(0,len(paths)):
			if (paths[iPath].vertices).shape[0] > maxPathLength:
				maxPathLength = (paths[iPath].vertices).shape[0]
				pathIndex = iPath
		
		frontPositions[fileToGet-1][iContour] =  paths[pathIndex].vertices




interpJetPositions   = np.zeros((numTimeSteps-1,len(frontSSHValues),nLon),dtype='float64')


for iT in range(0,numTimeSteps-1):
	for iJet in range(0,len(frontSSHValues)):
	
		interpJetPositions[iT,iJet,:] = griddata(frontPositions[iT][iJet][:,0], frontPositions[iT][iJet][:,1], lon, method='linear')

interpJetPositions = np.ma.masked_array(interpJetPositions,np.isnan(interpJetPositions))


fig = plt.figure(0)
ax  = fig.add_subplot(1,1,1)
ax.plot(sshValues[0:-1], aveAlongSSHGradientValues)


fig = plt.figure(1)
ax  = fig.add_subplot(1,1,1)
cPlot = ax.contourf(lon,lat,averagedSSH,sshValues,cmap=plt.cm.jet)
cbar = fig.colorbar(cPlot, ax=ax)


fig = plt.figure(3)
ax  = fig.add_subplot(1,1,1)
cPlot = ax.contourf(lon[5:-5],lat[5:-5],gradMag[5:-5,5:-5], 15, cmap=plt.cm.jet)
for iJet in range(0,len(frontSSHValues)):
	ax.plot(frontPositions[0][iJet][:,0],frontPositions[0][iJet][:,1],'k')
#for iPath in contourPoints:
#	plt.plot(iPath[:,0],iPath[:,1],'k')

cbar = fig.colorbar(cPlot, ax=ax)

fig = plt.figure(4)
ax  = fig.add_subplot(1,1,1)
cPlot = ax.contourf(lon[5:-5],lat[5:-5],D_SSH_Dlon[5:-5,5:-5], 15, cmap=plt.cm.jet)
cbar = fig.colorbar(cPlot, ax=ax)


contourID = 2
colourArray = ['k','b','m','g','r']
plt.figure(5)
plt.plot(frontPositions[0][contourID][:,0],frontPositions[0][contourID][:,1],colourArray[0])
plt.plot(lon,interpJetPositions[0,contourID,:],colourArray[4])


plt.show()
