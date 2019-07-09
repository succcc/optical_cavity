import numpy as np
import matplotlib.pyplot as plt
from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg
from sympy.physics.optics import *
from sympy import Symbol, Matrix
from cavity import Cavity
import sys

inf_meter = 10000000.0

khz = 1e3
mhz = 1e6
ghz = 1e9
thz = 1e12

c = 299792458.0

ms = 1 / khz
us = 1 / mhz
ns = 1 / ghz

m = 1
cm = 1e-2
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

p = BeamParameter(wavelen=1042*nm, z=0, w=c3.W0)
print("q(z) = ", p.q.n(), ", w(z) = ", p.w.n()/um)
p = FreeSpace(10*cm)*p
print("q(z) = ", p.q.n(), ", w(z) = ", p.w.n()/um)
p = ThinLens(10*cm)*p
p = FreeSpace(5*cm)*p
print("q(z) = ", p.q.n(), ", w(z) = ", p.w.n()/um)


class BeamPlot:
    def __init__(self):
        pg.setConfigOptions(antialias=False)
        self.app = QtGui.QApplication(sys.argv)
        self.win = pg.GraphicsWindow(title='Beam plotter')
        self.win.setWindowTitle('Beam plotter')
        # self.win.setGeometry(5, 115, 1910, 1070)

        bf_xlabels = [(0, '0'), (2048, '2048'), (4096, '4096')]
        bf_xaxis = pg.AxisItem(orientation='bottom')
        bf_xaxis.setTicks([bf_xlabels])

        bf_ylabels = [(0, '0'), (127, '128'), (255, '255')]
        bf_yaxis = pg.AxisItem(orientation='left')
        bf_yaxis.setTicks([bf_ylabels])

        self.beam = self.win.addPlot(
            title='Beam', axisItems={'bottom': bf_xaxis, 'left': bf_yaxis}
        )
        d1 = 20 * cm
        f1 = 20 * cm
        d2 = 40 * cm

        beam = BeamParameter(1042 * nm, 0, w=c3.W0)

        z_c = np.linspace(0, c3.d, 1000)


        z1 = np.linspace(0, d1, 1000)
        w = float(beam.w_0.n()) * np.sqrt(1 + (z1 / float(beam.z_r.n())) ** 2)
        self.beam.plot(-z1, w)
        self.beam.plot(-z1, -w)

        beam = ThinLens(f1)*FreeSpace(d1)*beam
        z2 = np.linspace(d1, d1+d2, 1000)
        w = float(beam.w_0.n()) * np.sqrt(1 + ((z2-d1+float(beam.z.n())) / float(beam.z_r.n())) ** 2)
        self.beam.plot(-z2, w)
        self.beam.plot(-z2, -w)
        print("q(z) = ", beam.q.n()/cm, ", w(z) = ", beam.w.n() / um)

        self.beam.showGrid(x=True, y=True)

    @staticmethod
    def start():
        if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
            QtGui.QApplication.instance().exec_()


if __name__ == '__main__':
    app = BeamPlot()
    app.start()
