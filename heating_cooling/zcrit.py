########################################################################
#
#
########################################################################

import numpy as np
from astropy import constants as const
from matplotlib import pyplot
pyplot.rcParams['ps.useafm'] = True
pyplot.rcParams['pdf.use14corefonts'] = True
pyplot.rcParams['text.usetex'] = True
#from matplotlib.font_manager import fontManager, FontProperties
#font=FontProperties(size='small')


Gamma = 1.0e-26  # erg/s
Z0 = 1.0 # relative to solar

Lambda0 = 1.2e-26 # cooling rate to T0, Z0
alpha = 1.0
T0 = 6000

sigma = 1.0 * 1e5 # 1 km/s
L = 20.0 * 3.08e18 # 10 pc
mu = 1.22 * const.u.cgs

def Pleq(T, Z):
    return Gamma * const.k_B.cgs * Z0 * T0**alpha * T**(1.0-alpha) / (Lambda0 * Z)

def PTeq(T, Z):
    return sigma**3 * mu * Z0 * (T0/T)**alpha / (L * Lambda0 * Z)


T = np.logspace(2.0, 4.2)
Z = 0.001

Pleq1 = Pleq(T, Z)

#pyplot.figure(figsize=(5,4))
pyplot.plot(T, Pleq1, color='k', ls='-')
pyplot.xscale('log')
pyplot.yscale('log')
pyplot.xlabel(r'Temperature (K)')
pyplot.ylabel(r'Pressure')
pyplot.axis([1e2,2e4,1e-12,1e-8])
#pyplot.legend(loc='lower left', prop=font)
pyplot.savefig('Pleq')
