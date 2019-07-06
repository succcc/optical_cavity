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

c1 = Cavity(roc1=-10*mm, roc2=-10*mm, pos1=0, d=10.001*mm, r1=0.99, r2=0.99, lamb=1042*nm)
c1.report()

c2 = Cavity(roc1=-100000000, roc2=-10*mm, pos1=0, d=5*mm, r1=0.99, r2=0.99, lamb=1042*nm)
c2.report()

p = BeamParameter(wavelen=1042*nm, z=100*mm, w=1*mm)
print(p.q)