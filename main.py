import numpy as np
import matplotlib.pyplot as plt
from sympy.physics.optics import *
from sympy import Symbol, Matrix
from cavity import Cavity

inf_meter = 10000000.0

khz = 1e3
mhz = 1e6
ghz = 1e9
thz = 1e12

c = 299792458.0

ms = 1 / khz
us = 1 / mhz
ns = 1 / ghz

mm = 1e-3
um = 1e-6
nm = 1e-9

ppm = 1e-6

# c1 = Cavity(roc1=-10*mm, roc2=-10*mm, pos1=0, d=10.001*mm, R1=0.98, R2=0.98, lamb=1042*nm, I0=1)
# c1Spectrum = c1.get_response(wCenter=1042 * nm, fWidth=30 * ghz, num=10000)
# plt.plot((c1Spectrum[0]-c/c1.lamb)/ghz,
#          c1Spectrum[1])
#
# c2 = Cavity(roc1=-inf_meter, roc2=-10.0*mm, pos1=0, d=5.0*mm, R1=0.9, R2=0.9, lamb=1042*nm, I0=1)
# c2Response = c2.get_response(wCenter=1042 * nm, fWidth=50 * ghz, num=10000)
# # plt.plot((c2Spectrum[0]-c/c2.lamb)/ghz,
# #          c2Spectrum[1])
# plt.plot((c2Response[0]-c/c2.lamb)/ghz, c2Response[2], label="Reflectance")
# plt.plot((c2Response[0]-c/c2.lamb)/ghz, c2Response[3], label="Transmission")
#
# plt.xlabel("Frequency, zero at 1042 nm (GHz)")
# plt.ylabel("I/I0")
#
# plt.legend()
# plt.tight_layout()
# plt.show()


L = 10*ppm
T = 90*ppm
R = 1-L-T
c3 = Cavity(roc1=-inf_meter, roc2=-10*mm, pos1=0, d=5.0*mm, R1=R, R2=R, L1=L, L2=L, lamb=1042*nm, I0=1)
c3.report()
