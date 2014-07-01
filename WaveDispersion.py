import numpy as np
import os
from scipy.io import netcdf as scipy_netcdf
from matplotlib import pyplot as plt

def GetMeridionalWaveNumber(zonalWaveNumber,U_Profile,yIndex,beta,deltaY):
	
	#Compute the curvature of the zonal flow profile
	deltaY_m = deltaY*1.0e3
	d2Udy2 = (U_Profile[2::] + U_Profile[0:-2] - 2.0*U_Profile[1:-1]) / (deltaY_m*deltaY_m)
	
	modBeta = beta - d2Udy2[yIndex]
	meridionalWaveNumber = np.sqrt((modBeta/U_Profile[yIndex]) - (zonalWaveNumber*zonalWaveNumber)  )
	
	return meridionalWaveNumber	
	


meanDepth = 4000.0
f0   = -1.0e-4
beta = 1.5e-11

outputDirectory = '/home/chris/GOLD/IsoTopo_1_Processing/'

modelOutputFileName = 'postprocessed_IsoTopo1.nc'
topoFileName   = 'isoTopo.nc'

modelFile = scipy_netcdf.netcdf_file(outputDirectory + modelOutputFileName, 'r')
topoFile  = scipy_netcdf.netcdf_file(outputDirectory + topoFileName, 'r')
#Get dimensions
xName = 'x'
yName = 'y'
zlName = 'zl'
ziName = 'zi'

x  = modelFile.variables[xName][:]
y  = modelFile.variables[yName][:]
zl = modelFile.variables[zlName][:]
zi = modelFile.variables[ziName][:]

deltaX = x[1]-x[0]
deltaY = y[1]-y[0]

uName    = 'u_time_mean'
vName    = 'v_time_mean'
ekeName  = 'eke'
topoName = 'topo'

uTimeMean   = modelFile.variables[uName][:,:,:]
vTimeMean   = modelFile.variables[vName][:,:,:]
ekeTimeMean = modelFile.variables[ekeName][:,:,:]
topo        = meanDepth - topoFile.variables[topoName][:,:]
#zonalWaveNumber,U_Profile,yIndex,beta,deltaY

GetMeridionalWaveNumber(9,uTimeMean[0,:,x.shape[0]/2+20],y.shape[0]/3,beta,deltaY)


modelFile.close()
topoFile.close()

