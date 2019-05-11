# -*- coding: utf-8 -*-
# Program to print out the energies of elementary excitations
# of transvers field ising model.

import math
import numpy as np
from exact import H_TFIM, get_lowest_modes, get_options

def main():
    # read the parameters from input
    # L: the amount of sites
    # J, g: parameters of H_{TFIM}
    L, J, g = get_options()

    # get the discrete k-values; eq. (2.67)
    n = np.array(range(L))
    k_even = (2 * math.pi / L) * n
    k_odd = (math.pi / L) * (2 * n + 1)

    # get the energies of the elementary excitations; eq. (2.74)
    k = np.concatenate((k_even, k_odd))
    E = - 2 * J * np.sqrt(g**2 + 1 - 2 * g * np.cos(k))
    print('E(k):', np.sort(E))
    # Remark that for g=0 elementary excitations are created in pairs,
    # while for g>>J single elementary excitations are created.
    # This can be observed by comparison with exact diagonalization.
    # e.g. L=8, J=2 and g=0 VS g=50
    


        

# don't run main if included 
if __name__ == "__main__":
    main()