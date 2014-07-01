import numpy as np
import os
from scipy.io import netcdf as scipy_netcdf
from matplotlib import pyplot as plt
from scipy.signal import firwin, filtfilt, lfilter
from scipy.fftpack import fft, fftfreq, fftshift, ifft, fft2, ifft2


topoLengthByExp = np.asarray([50,50,100,50,100,25,200,150,100,100])
topoHeigthByExp = np.asarray([1000,2000,2000.,2000,2000,2000,2000,2000,1000,3000])
topoWidthByExp  = np.asarray([100,100,100,50,50,50,50,50,50,50])

meanderWavelength = np.asarray([667.73339303,516.80004621,573.86671797,540.80004835,561.06671683,
					 699.73339589,550.40004921,486.40004349,580.26671855,573.86671797])
					 
intEKE  = np.asarray([9.53001,10.5346,9.28285,9.64423,10.4709,11.0303,10.1198,10.302,9.89645,10.7778])
stormTrackLength = np.asarray([1638.40014648, 1939.20017338, 1670.40014935, 1926.40017223,
					1606.40014362,2227.20019913,1523.20013618,1574.40014076, 1689.60015106, 1689.60015106])

vAmp = np.asarray([0.96,0.8,0.85,0.81,0.6,0.69,0.4,0.43,0.5,0.55])
peakU = np.asarray([0.32,0.32,0.3,0.36,0.32,0.32,0.32,0.3,0.32,0.32])

meanDepth = 4000.0
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

			
sortIndex = np.argsort(topoLengthByExp[3:8])

print sortIndex
print topoLengthByExp[3:8][sortIndex]
print type(sortIndex)
print topoHeigthByExp[8],topoHeigthByExp[5],topoHeigthByExp[9]
fig = plt.figure(1)
ax  = fig.add_subplot(1,1,1)
ax.plot(topoLengthByExp[3:8][sortIndex],stormTrackLength[3:8][sortIndex],'k*')
ax.plot(topoLengthByExp[3:8][sortIndex],stormTrackLength[3:8][sortIndex],'k-')
ax.set_ylabel('Storm Track Length (km)')
ax.set_xlabel('Zonal Topo Length (km)')



print [topoHeigthByExp[5],topoHeigthByExp[8:10]]
fig = plt.figure(2)
ax  = fig.add_subplot(1,1,1)
ax.plot([topoHeigthByExp[8],topoHeigthByExp[5],topoHeigthByExp[9]],[intEKE[8],intEKE[5],intEKE[9]],'k*')
ax.plot([topoHeigthByExp[8],topoHeigthByExp[5],topoHeigthByExp[9]],[intEKE[8],intEKE[5],intEKE[9]],'k-')
ax.set_ylabel('$\int$EKE (m$^{4}$s$^{-2}$)')
ax.set_xlabel('Topo Height (m)')


fig = plt.figure(3)
ax  = fig.add_subplot(1,1,1)
ax.plot([topoHeigthByExp[8],topoHeigthByExp[5],topoHeigthByExp[9]],[stormTrackLength[8],stormTrackLength[5],stormTrackLength[9]],'k*')
ax.plot([topoHeigthByExp[8],topoHeigthByExp[5],topoHeigthByExp[9]],[stormTrackLength[8],stormTrackLength[5],stormTrackLength[9]],'k-')
ax.set_ylabel('Storm Track Length (km)')
ax.set_xlabel('Topo Height (m)')

fig = plt.figure(4)
ax  = fig.add_subplot(1,1,1)
ax.plot([topoHeigthByExp[8],topoHeigthByExp[5],topoHeigthByExp[9]],[vAmp[8],vAmp[5],vAmp[9]],'k*')
ax.plot([topoHeigthByExp[8],topoHeigthByExp[5],topoHeigthByExp[9]],[vAmp[8],vAmp[5],vAmp[9]],'k-')
ax.set_ylabel('Standing Wave Amplitude (m/s)')
ax.set_xlabel('Topo Height (m)')


fig = plt.figure(5)
ax  = fig.add_subplot(2,1,1)
ax.annotate('a',xy=(0.01,1.025),xycoords='axes fraction',fontsize=15)

ln1= ax.plot(topoLengthByExp[3:8][sortIndex]*6.4,vAmp[3:8][sortIndex],'k*',label='Wave Amplitude')
ax.plot(topoLengthByExp[3:8][sortIndex]*6.4,vAmp[3:8][sortIndex],'k-')
ax2 = ax.twinx()
ln2=ax2.plot(topoLengthByExp[3:8][sortIndex]*6.4,stormTrackLength[3:8][sortIndex],'ko',label='Storm Track Length')
ax2.plot(topoLengthByExp[3:8][sortIndex]*6.4,stormTrackLength[3:8][sortIndex],'k-')

ax.set_ylabel('Standing Wave Amplitude (m/s)')
ax2.set_ylabel('Storm Track Length (km)')
ax.get_xaxis().set_visible(False)

lns = ln1+ln2
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc='upper right', shadow=True, prop={'size':10})

ax1  = fig.add_subplot(2,1,2)
ax1.annotate('b',xy=(0.01,1.025),xycoords='axes fraction',fontsize=15)

ln3=ax1.plot(topoLengthByExp[3:8][sortIndex]*6.4,(peakU[3:8][sortIndex]),'k^',label='Max U')
ax1.plot(topoLengthByExp[3:8][sortIndex]*6.4,(peakU[3:8][sortIndex]),'k-')
ax2 = ax1.twinx()
ln4=ax2.plot(topoLengthByExp[3:8][sortIndex]*6.4,densityLayers[1]*intEKE[3:8][sortIndex]/1000.0,'k+',label='Averaged EKE')
ax2.plot(topoLengthByExp[3:8][sortIndex]*6.4,densityLayers[1]*intEKE[3:8][sortIndex]/1000.0,'k-')
lns = ln3+ln4

ax1.set_ylabel('Maximum U (m/s)')
ax2.set_ylabel('Integrated EKE (kJ.m$^{-1}$)')

labs = [l.get_label() for l in lns]
ax1.legend(lns, labs, loc='upper right', shadow=True, prop={'size':10})

ax1.set_xlabel('Zonal Topo Length (km)')

plt.show()					
					




