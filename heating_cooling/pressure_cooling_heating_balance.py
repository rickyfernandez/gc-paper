########################################################################
#
# Cooling rate example script
#
#
# Copyright (c) 2013, Enzo/Grackle Development Team.
#
# Distributed under the terms of the Enzo Public Licence.
#
# The full license is in the file LICENSE, distributed with this 
# software.
########################################################################

import copy
import numpy as np
import sys

import h5py

from pygrackle.grackle_wrapper import *
from pygrackle.fluid_container import FluidContainer

from utilities.api import \
     setup_fluid_container, \
     set_cosmology_units, \
     get_cooling_units

from utilities.physical_constants import \
     gravitational_constant_cgs,\
     mass_hydrogen_cgs, \
     sec_per_Gyr, \
     cm_per_mpc, \
     boltzmann_constant_cgs

def cooling_to_freefall_map(dlo, dhi, tlo, thi, num_points, my_chemistry, z_frac):

    # Call convenience function for setting up a fluid container.
    # This container holds the solver parameters, units, and fields.
    temperature  = np.logspace(tlo, thi, num_points)
    density      = np.logspace(dlo, dhi, num_points)*mass_hydrogen_cgs
    ratio_values = np.zeros((density.size, density.size))
    Temp  = np.zeros(num_points)

    for i, d in enumerate(density):

        fc = setup_fluid_container(my_chemistry,
                                   density=d,
                                   metal_mass_fraction=z_frac*0.02041,
                                   temperature=temperature,
                                   current_redshift=0,
                                   converge=False)

        calculate_temperature(fc, 1.0)
        calculate_cooling_time(fc, 1.0)
        t_sort = np.argsort(fc["temperature"])

        temp = fc["temperature"][t_sort]
        Temp[i] = temp[ np.argmax(fc["cooling_time"][t_sort]) ]

    return density, Temp


if __name__ == "__main__":

    # Set solver parameters
    my_chemistry = chemistry_data()
    my_chemistry.use_grackle = 1
    my_chemistry.with_radiative_cooling = 0
    my_chemistry.primordial_chemistry = 0
    my_chemistry.metal_cooling = 1
    my_chemistry.UVbackground = 1
    my_chemistry.photoelectric_heating = 1
    my_chemistry.grackle_data_file = "CloudyData_UVB=HM2012.h5"
    #my_chemistry.grackle_data_file = "CloudyData_noUVB.h5"

    # Set units
    my_chemistry.comoving_coordinates = 0 # proper units
    my_chemistry.a_units = 1.0
    my_chemistry.density_units = mass_hydrogen_cgs # rho = 1.0 is 1.67e-24 g
    my_chemistry.length_units = cm_per_mpc         # 1 Mpc in cm
    my_chemistry.time_units = sec_per_Gyr          # 1 Gyr in s
    my_chemistry.velocity_units = my_chemistry.length_units / my_chemistry.time_units


    # density range in log space
    dlo = 1.0
    dhi = 9.0

    # temperature range in log space
    tlo = 1.0
    thi = 9.0

    # fraction of solar metallicity
    z = 0.01

    # resolution - number of bins in density/temperature
    # dimension
    num_points = 1024

    # create equlibrium curve
    dens, temp = cooling_to_freefall_map(dlo, dhi, tlo, thi, num_points, my_chemistry, z)

    # store the equilibrum curve
    f = h5py.File("equilibrium_curve.hdf5", "w")
    f['/denstity'] = dens
    f['/temperature'] = temp
    f.close()
