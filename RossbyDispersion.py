import numpy as np
from matplotlib import pyplot as plt
from scipy.fftpack import fft, fftfreq, fftshift, ifft, fft2, ifft2

nX = 800
deltaX = 6.4e3
U = 1.0 
BETA = 1.5e-11
Rd = 25.0e3
Kd = 1.0/Rd
wavelengthY = 500.0e3

BETA_KM_DAY = BETA * (60.0*60.0*24.0) / (1000.0)


x = np.arange(0,deltaX*nX,deltaX)

kx = fftshift(fftfreq(nX,d=deltaX))
ky = 1.0/wavelengthY

cp = U - ( BETA/((kx*kx + ky*ky) + (Kd*Kd)) )


fig = plt.figure(1)
ax  = fig.add_subplot(1,1,1)
ax.plot((1.0/kx[nX/2::])*1.0e-3, cp[nX/2::]  )
plt.show()

