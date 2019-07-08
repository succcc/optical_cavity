import numpy as np

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


class Cavity:
    def __init__(self, roc1, roc2, pos1, d, r1, r2, lamb, i0):
        self.pos1 = pos1
        self.pos2 = pos1 + d
        self.roc1 = roc1
        self.roc2 = roc2
        self.r1 = r1
        self.r2 = r2
        self.d = d
        self.lamb = lamb
        self.i0 = i0

        self.nu0 = c / self.lamb
        self.W0 = np.power((self.lamb / np.pi) ** 2 * (
                -self.d * (self.roc1 + d) * (self.roc2 + self.d) * (self.roc1 + self.roc2 + self.d))
                           / (self.roc1 + self.roc2 + 2 * self.d) ** 2, 1 / 4)
        self.z0 = np.pi * self.W0 ** 2 / self.lamb
        self.FSR = 3e8 / 2 / self.d
        self.g1g2 = (1 + self.d / self.roc1) * (1 + self.d / self.roc2)
        self.finesse = np.pi * (r1 * r2) ** (1 / 4) / (1 - np.sqrt(r1 * r2))
        self.a_s = 0  # Assuming no attenuation in cavity, only consider reflectivity of mirrors

        self.ar = self.a_s + 1 / (2 * d) * np.log(1 / (r1 * r2))
        self.absr = np.sqrt(np.exp(-2*self.ar*self.d))

        self.finesse2 = np.pi * np.exp(-self.ar * d / 2) / (1 - np.exp(-self.ar * d))
        self.delta_nu = self.FSR / self.finesse  # approx
        self.tau_p = 1 / (2 * np.pi * self.delta_nu)
        self.tau_p2 = 1 / (c * self.ar)
        self.qFactor = self.nu0 / self.delta_nu
        self.qFactor_correct = 2 * np.pi * self.nu0 / (c * self.ar)

        self.imax = self.i0/((1-self.absr)**2)

    def intensity(self, freq):
        return self.imax/(1 + (2*self.finesse/np.pi)**2 * np.sin(np.pi*freq/self.FSR)**2)

    def setD(self, d):
        self.__init__(self.roc1, self.roc2, self.pos1, d, self.r1, self.r2, self.lamb)

    def setPos1(self, pos1):
        self.__init__(self.roc1, self.roc2, pos1, self.d, self.r1, self.r2, self.lamb)

    def getSpectrum(self, wStart=None, wEnd=None, fStart=None, fEnd=None, wCenter=None, fWidth=None, num=1000):
        isWavelength = wStart or wEnd
        isFrequencyIntv = fStart or fEnd
        isWavelengthCenter = wCenter or fWidth
        if isWavelength and not (isWavelengthCenter or isFrequencyIntv):
            lamb = np.linspace(wStart, wEnd, num)
            freq = c / lamb
        elif isFrequencyIntv and not (isWavelength or isWavelengthCenter):
            freq = np.linspace(fStart, fEnd, num)
        elif isWavelengthCenter and not (isWavelength or isFrequencyIntv):
            freq = np.linspace(c/wCenter - fWidth, c/wCenter + fWidth, num)
        else:
            print("Please check the input variables...")
            return 0, 0
        return freq, self.intensity(freq)

    def report(self):
        print(c)
        print("======================================")
        print("Printing parameters and derived values")
        print("-------------Resonance----------------")

        print("Q = {:0.2f}".format(self.qFactor))
        print("F = {:0.2f}".format(self.finesse))
        print("FSR = {:0.2f} GHz".format(self.FSR / ghz))
        print("FWHM = {:0.2f} MHz".format(self.delta_nu / mhz))
        print("nu0 = {:0.2f} THz".format(self.nu0 / thz))

        print("--------------Gaussian----------------")
        print("W0 = {:0.2f} um".format(self.W0 / um))
        print("z0 = {:0.2f} mm".format(self.z0 / mm))
        print("g1g2 = {:0.2f}".format(self.g1g2))

        print("")