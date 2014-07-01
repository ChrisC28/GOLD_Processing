import numpy as np
from matplotlib import pyplot as plt

def ComputeMontPot(upperInterfaceHeight,thickness,reducedGravs, densities, rho0,g0):
	
	nZ_l = thickness.shape[1]
	
	
	montPot = np.zeros(thickness.shape,dtype='float64')
	
	for iLayer in range(0,nZ_l):
		
		if iLayer == 0:
		
			montPot[:,iLayer,:,:] = (densities[iLayer] * g0 * upperInterfaceHeight)/rho0
		
		else:
			
			montPot[:,iLayer,:,:] = (montPot[:,iLayer-1,:,:] + \
			                        ( reducedGravs[iLayer] * thickness[:,iLayer,:,:] ) ) 
			
	return montPot 
	
def ComputeEPvector(uEddy,vEddy,montPot, uMean, vMean, refThickness, densityLayers, deltaX, deltaY, g0,f0,rho0):
	
	nX      = uEddy.shape[3]
	nY      = uEddy.shape[2]
	nLayers = uEddy.shape[1]
	nT      = uEddy.shape[0]
	
	EP_X = np.zeros((nT,nLayers,nX-2),dtype='float64')
	EP_Y = np.zeros((nT,nLayers,nX-2),dtype='float64')
	EP_Z = np.zeros((nT,nLayers,nX-2),dtype='float64')
	
	
	
	#print nLayers
	#print refThickness
	
	uEddyT   = np.squeeze(uEddy[:,:,1:nY-1,1:nX-1])
	vEddyT   = np.squeeze(vEddy[:,:,1:nY-1,1:nX-1])
	#montPotT = montPot[:,:,1:nY-1,1:nX-1]
	psi = np.squeeze(montPot[:,:,1:nY-1,1:nX-1])/f0
	
	
	dudx = np.squeeze((uEddy[:,:,1:nY-1,2:nX] - uEddy[:,:,1:nY-1,0:nX-2]) / (2.0*deltaX))
	dvdx = np.squeeze((vEddy[:,:,1:nY-1,2:nX] - vEddy[:,:,1:nY-1,0:nX-2]) / (2.0*deltaX))

	dudy = np.squeeze((uEddy[:,:,2:nY,1:nX-1] - uEddy[:,:,0:nY-2,1:nX-1]) / (2.0*deltaY))
	dvdy = np.squeeze((vEddy[:,:,2:nY,1:nX-1] - vEddy[:,:,0:nY-2,1:nX-1]) / (2.0*deltaY))

	dudp = np.squeeze(vertDiff(uEddyT,densityLayers))
	dvdp = np.squeeze(vertDiff(vEddyT,densityLayers))

	dPsidp = vertDiff(psi,densityLayers)
	
	EP_X[:,:,:] = (uMean[:,:,1:uMean.shape[-1]-1] * ( vEddyT * vEddyT - psi * dvdx)) + (vMean[:,:,1:vMean.shape[-1]-1] * (-uEddyT * vEddyT + psi * dudx))
	EP_Y[:,:,:] = (uMean[:,:,1:uMean.shape[-1]-1] * (-uEddyT * vEddyT - psi * dvdy)) + (vMean[:,:,1:vMean.shape[-1]-1] * ( uEddyT * uEddyT + psi * dudy))
	EP_Z[:,:,:] = (uMean[:,:,1:uMean.shape[-1]-1] * (vEddyT * dPsidp -  psi * dvdp)) + (vMean[:,:,1:vMean.shape[-1]-1] * (-uEddyT * dPsidp + psi * dudp)) 
	
	
	for iLayer in range(0,nLayers):
		EP_Z[:,iLayer,:] = f0*densityLayers[iLayer]/(g0*refThickness[iLayer]) * EP_Z[:,iLayer,:]
		
	return EP_X, EP_Y, EP_Z

def ComputeEP_Components(uEddy,vEddy,montPot, uMean, vMean, refThickness, densityLayers, deltaX, deltaY, g0,f0,rho0):
	
	nX      = uEddy.shape[3]
	nY      = uEddy.shape[2]
	nLayers = uEddy.shape[1]
	nT      = uEddy.shape[0]
	
	EP_X = np.zeros((nT,nLayers,nX-2),dtype='float64')
	EP_Y = np.zeros((nT,nLayers,nX-2),dtype='float64')
	EP_Z = np.zeros((nT,nLayers,nX-2),dtype='float64')
	
	
	
	#print nLayers
	#print refThickness
	
	uEddyT   = np.squeeze(uEddy[:,:,1:nY-1,1:nX-1])
	vEddyT   = np.squeeze(vEddy[:,:,1:nY-1,1:nX-1])
	#montPotT = montPot[:,:,1:nY-1,1:nX-1]
	psi = np.squeeze(montPot[:,:,1:nY-1,1:nX-1])/f0
	
	
	dudx = np.squeeze((uEddy[:,:,1:nY-1,2:nX] - uEddy[:,:,1:nY-1,0:nX-2]) / (2.0*deltaX))
	dvdx = np.squeeze((vEddy[:,:,1:nY-1,2:nX] - vEddy[:,:,1:nY-1,0:nX-2]) / (2.0*deltaX))

	dudy = np.squeeze((uEddy[:,:,2:nY,1:nX-1] - uEddy[:,:,0:nY-2,1:nX-1]) / (2.0*deltaY))
	dvdy = np.squeeze((vEddy[:,:,2:nY,1:nX-1] - vEddy[:,:,0:nY-2,1:nX-1]) / (2.0*deltaY))

	dudp = np.squeeze(vertDiff(uEddyT,densityLayers))
	dvdp = np.squeeze(vertDiff(vEddyT,densityLayers))

	dPsidp = vertDiff(psi,densityLayers)
	
	reynoldsStress = (uMean[:,:,1:uMean.shape[-1]-1] * ( vEddyT * vEddyT)) + (vMean[:,:,1:vMean.shape[-1]-1] * (-uEddyT * vEddyT))
	potentialFlux  = (uMean[:,:,1:uMean.shape[-1]-1] * ( - psi * dvdx)) + (vMean[:,:,1:vMean.shape[-1]-1] * (psi * dudx))

	return reynoldsStress, potentialFlux
	#EP_Y[:,:,:] = (uMean[:,:,1:uMean.shape[-1]-1] * (-uEddyT * vEddyT - psi * dvdy)) + (vMean[:,:,1:vMean.shape[-1]-1] * ( uEddyT * uEddyT + psi * dudy))

	
	
	
def vertDiff(inputField, densities):
	nT, nLayers, nX = inputField.shape
	diffField = np.zeros(inputField.shape,dtype='float64')
	
	for iLayer in range(0,nLayers):
		
		if 0 == iLayer:
		
			diffField[:,iLayer,:] = (inputField[:,iLayer,:] - inputField[:,iLayer+1,:]) / (densities[iLayer] - densities[iLayer+1]) 
		
		elif nLayers-1 == iLayer:
		
			diffField[:,iLayer,:] = (inputField[:,iLayer-1,:] - inputField[:,iLayer,:]) / (densities[iLayer-1] - densities[iLayer])
		
		else: 
		
			diffField[:,iLayer,:] = 0.5*(inputField[:,iLayer+1,:] - inputField[:,iLayer-1,:]) / (densities[iLayer+1] - densities[iLayer-1])
			
	return diffField		
			
			
	
	
	 
	 
	
	
	
	
	
	
	

