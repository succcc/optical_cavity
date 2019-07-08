import numpy as np
import matplotlib.pyplot as plt
from sympy.physics.optics import *
from sympy import Symbol, Matrix
from cavity import Cavity
import cavity


khz = 1e3
mhz = 1e6
ghz = 1e9
thz = 1e12

c = 299792458

ms = 1 / khz
us = 1 / mhz
ns = 1 / ghz

mm = 1e-3
um = 1e-6
nm = 1e-9

c1 = Cavity(roc1=-10*mm, roc2=-10*mm, pos1=0, d=10.001*mm, r1=0.95, r2=0.95, lamb=1042*nm, i0=1)
c1Spectrum = c1.getSpectrum(wCenter=1042 * nm, fWidth=30 * ghz, num=100000)
plt.plot((c1Spectrum[0]-c/c1.lamb)/ghz,
         c1Spectrum[1])

c2 = Cavity(roc1=-10000000.0, roc2=-10.0*mm, pos1=0, d=5.0*mm, r1=0.99, r2=0.99, lamb=1042*nm, i0=1)
c2Spectrum = c2.getSpectrum(wCenter=1042 * nm, fWidth=30 * ghz, num=100000)
plt.plot((c2Spectrum[0]-c/c2.lamb)/ghz,
         c2Spectrum[1])

plt.xlabel("Frequency, zero at 1042 nm (GHz)")
plt.ylabel("Intensity, same input")
plt.show()
