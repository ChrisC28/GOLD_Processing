import numpy as np
import os
from scipy.io import netcdf as scipy_netcdf
from matplotlib import pyplot as plt


baseDir = '/home/chris/GOLD/IsoTopo4/'
outputDirBaseString = 'output'
startYear = 50
endYear   = 51



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


	
fileObject = scipy_netcdf.netcdf_file(uvFilesList[0], 'r')
u  = fileObject.variables['u'][:,:,:,:]
v  = fileObject.variables['v'][:,:,:,:]
e  = fileObject.variables['e'][:,:,:,:]
	
iTime = 0	
iLayer_uv = 0
iLayer_e  = 4

fig = plt.figure(1)
ax  = fig.add_subplot(1,1,1)
cs = ax.contourf(xq,yh,u[iTime,iLayer_uv,:,:],15,cmap=plt.cm.RdBu_r)
fig.colorbar(cs,ax=ax)

fig = plt.figure(2)
ax  = fig.add_subplot(1,1,1)
cs = ax.contourf(xh,yq,v[iTime,iLayer_uv,:,:],15,cmap=plt.cm.RdBu_r)
fig.colorbar(cs,ax=ax)


fig = plt.figure(3)
ax  = fig.add_subplot(1,3,1,aspect='equal')
cs = ax.contourf(xh[nX_h/2::],yh,e[iTime,iLayer_e,:,nX_h/2::],30,cmap=plt.cm.jet)
ax  = fig.add_subplot(1,3,2,aspect='equal')
cs = ax.contourf(xh[nX_h/2::],yh,e[iTime+2,iLayer_e,:,nX_h/2::],30,cmap=plt.cm.jet)
ax  = fig.add_subplot(1,3,3,aspect='equal')
cs = ax.contourf(xh[nX_h/2::],yh,e[iTime+5,iLayer_e,:,nX_h/2::],30,cmap=plt.cm.jet)
plt.show()


 
	
	
	
fileObject.close()

