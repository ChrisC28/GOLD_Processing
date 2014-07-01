import numpy as np
import os
from scipy.io import netcdf as scipy_netcdf
from matplotlib import pyplot as plt
from matplotlib.path import Path
from scipy.signal import firwin, filtfilt, lfilter, find_peaks_cwt
from scipy.fftpack import fft, fftfreq, fftshift, ifft, fft2, ifft2


xs = np.arange(0, np.pi, 0.05)
data = np.sin(xs)


peakind = find_peaks_cwt(data, np.arange(1,5))

print xs[peakind], data[peakind]
plt.figure(1)
plt.plot(xs,data)
plt.plot(xs[peakind], data[peakind],'b*')
plt.show()
