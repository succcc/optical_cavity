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
    def __init__(self, roc1, roc2, pos1, d, R1, R2, lamb, I0):
        self.pos1 = pos1
        self.pos2 = pos1 + d
        self.roc1 = roc1
        self.roc2 = roc2
        self.R1 = R1
        self.R2 = R2
        self.T1 = 1-R1
        self.T2 = 1-R2
        self.d = d
        self.lamb = lamb
        self.I0 = I0
        self.I0_in = self.T1*I0  # Initial intensity inside cavity
        self.E0 = np.sqrt(I0)

        self.nu0 = c / self.lamb
        self.W0 = np.power((self.lamb / np.pi) ** 2 * (
                -self.d * (self.roc1 + d) * (self.roc2 + self.d) * (self.roc1 + self.roc2 + self.d))
                           / (self.roc1 + self.roc2 + 2 * self.d) ** 2, 1 / 4)
        self.z0 = np.pi * self.W0 ** 2 / self.lamb
        self.FSR = c / 2 / self.d
        self.g1g2 = (1 + self.d / self.roc1) * (1 + self.d / self.roc2)
        self.finesse = np.pi * (R1 * R2) ** (1 / 4) / (1 - np.sqrt(R1 * R2))
        self.a_s = 0  # Assuming no attenuation in cavity, only consider reflectivity of mirrors

        self.ar = self.a_s + 1 / (2 * d) * np.log(1 / (R1 * R2))
        self.r_abs = np.sqrt(np.exp(-2*self.ar*self.d))

        self.finesse2 = np.pi * np.exp(-self.ar * d / 2) / (1 - np.exp(-self.ar * d))
        self.delta_nu = self.FSR / self.finesse  # approx
        self.tau_p = 1 / (2 * np.pi * self.delta_nu)
        self.tau_p2 = 1 / (c * self.ar)
        self.qFactor = self.nu0 / self.delta_nu
        self.qFactor_correct = 2 * np.pi * self.nu0 / (c * self.ar)

        self.I_max = self.I0_in/((1-self.r_abs)**2)

    def intensity_i(self, freq):
        return self.I_max/(1 + (2*self.finesse/np.pi)**2 * np.sin(np.pi*freq/self.FSR)**2)

    def set_d(self, d):
        self.__init__(self.roc1, self.roc2, self.pos1, d, self.R1, self.R2, self.lamb, self.I0)

    def set_pos1(self, pos1):
        self.__init__(self.roc1, self.roc2, pos1, self.d, self.R1, self.R2, self.lamb, self.I0)

    def get_response(self, wStart=None, wEnd=None, fStart=None, fEnd=None, wCenter=None, fWidth=None, num=1000):
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
        return freq, self.intensity_i(freq)

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