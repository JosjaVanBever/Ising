# -*- coding: utf-8 -*-
# Program to analyse the gap of transvers field ising model.

from exact import H_TFIM, get_lowest_modes, get_order_parameter
import numpy as np

def main():
    # set printing format for numpy
    np.set_printoptions(precision=4)

    # parameters in the hamiltonian
    J = 1.0
    g = 0.5

    # keep track of the energy gaps
    gaps = []  # gap
    Ls = []    # corresponding L
    # scan over lattice size L
    for L in range(2,13):
        Ls.append(L)
        ham = H_TFIM(L,J,g)
        val,vec = get_lowest_modes(ham, 2)
        gaps.append(val[1]-val[0])

    # the gap closes towards large L
    # this demonstrates the statement in section 2.1.2
    print("gap between lowest 2 modes as function of L:")
    for i in range(len(gaps)):
        print('%d: %.5f'% (Ls[i], gaps[i]))

    print()
    # fix now the lattice size L
    L = 8
    J = 20.

    # keep track of the gap to excitations
    gaps = []  # gap
    gs = []    # corresponding g
    ms = []    # corresponding m
    # scan over parameter g of the hamiltonian
    for g in np.linspace(0, 2, 15):
        gs.append(g)
        ham = H_TFIM(L,J,g)
        val,vec = get_lowest_modes(ham, 3)
        # get the gap and m
        ms.append(get_order_parameter(vec[:,0], ham))
        # print('m = {0:.3f}'.format(get_order_parameter(vec[:,0], ham)))
        gaps.append(val[2]-val[1])

    print("gap to excitation and m as function of g:")
    for i in range(len(gaps)):
        print('g=%.2f: E=%.5f, m=%.3f'% (gs[i], gaps[i], ms[i]))


# don't run main if included 
if __name__ == "__main__":
    main()
