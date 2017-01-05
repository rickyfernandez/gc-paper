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

chi = 1.0

sigma = 2.5 * 1e5 # 1 km/s
L = 5.0 * 3.08e18 # 10 pc
mu = 1.22 * const.m_p.cgs
c_s = 0.7*10.0 *(6000/10000.0)**0.5 * 1e5 # 10 km/s
MBE = 1.0e6 * 2e33 # 10^6 msun

def Pleq(T, Z):
    return Gamma * const.k_B.cgs * Z0 * T0**alpha * T**(1.0-alpha) / (Lambda0 * Z)

def PTeq(T, Z):
    return sigma**3 * mu * const.k_B.cgs * Z0 * T0**alpha * T**(1.0-alpha) / (L * Lambda0 * Z)

def PBE():
    return 19.9 * chi * c_s**8 / (const.G.cgs**3 * MBE**2)

def ZcritL(Gamma, MBE, T):
    return Z0 * const.k_B.cgs * T0**alpha * const.G.cgs**3 * T**(1.0-alpha) * Gamma * MBE**2 / (19.9 * Lambda0 * chi * c_s**8)

def ZcritT(MBE, T):
    return Z0 * sigma**3 * mu * const.k_B.cgs * T0**alpha * const.G.cgs**3 * T**(1.0-alpha) * MBE**2 / (19.9 * Lambda0 * L * chi * c_s**8)
     
T = np.logspace(2.0, 4.2)
Z = 0.001

Pleq1 = Pleq(T, Z)
PTeq1 = PTeq(T, Z)

print PBE()

pyplot.figure(figsize=(5,4))
pyplot.plot(T, Pleq1, color='k', ls='-', label='PLeq')
pyplot.plot(T, PTeq1, color='k', ls='--', label='PTeq')
pyplot.xscale('log')
pyplot.yscale('log')
pyplot.xlabel(r'Temperature (K)')
pyplot.ylabel(r'Pressure')
pyplot.axis([1e2,2e4,1e-12,1e-8])
#pyplot.legend(loc='lower left', prop=font)
pyplot.savefig('Peq')


Gamma = np.logspace(-28, -24)

ZcritL1 = ZcritL(Gamma, MBE, T0)
ZcritT1 = np.full_like(ZcritL1, ZcritT(MBE, T0))

pyplot.figure(figsize=(5,4))
pyplot.plot(Gamma, ZcritL1, color='k', ls='-', label='ZcritL')
pyplot.plot(Gamma, ZcritT1, color='k', ls='--', label='ZcritT')
pyplot.xscale('log')
pyplot.yscale('log')
pyplot.xlabel(r'$\Gamma$ (erg s cm$^{-3}$)')
pyplot.ylabel(r'$Z_{\rm crit}$')
pyplot.axis([1e-28,1e-24,1e-4,1])
#pyplot.legend(loc='lower left', prop=font)
pyplot.savefig('Zcrit')


