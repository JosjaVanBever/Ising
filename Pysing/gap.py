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
    gaps1 = [] # gap between "ground and 1st excitation"
    gaps2 = [] # gap between "1st and 2nd excitation"
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
        gaps1.append(val[1]-val[0])
        gaps2.append(val[2]-val[1])

    print("gap to excitations and m as function of g:")
    for i in range(len(gaps)):
        print('g=%.2f: E(1-0)=%.5f, E(2-1)=%.5f, m=%.3f'% \
            (gs[i], gaps1[i], gaps2[i], ms[i]))

    # m as function of L for g around g_c
    gs = []
    ms = []
    Ls = range(2,13)
    for g in np.linspace(0.8, 1.2, 8):
        gs.append(g)
        Lms = []
        for L in Ls:
            ham = H_TFIM(L,J,g)
            val,vec = get_lowest_modes(ham, 1)
            Lms.append(get_order_parameter(vec[:,0], ham))
        ms.append(Lms)

    for i in range(len(gs)):
        print('g=%.2f:'% (gs[i]), end = ' ')
        for j in range(len(ms[i])):
            print('(L=%d: m=%.3f)'% (Ls[j], ms[i][j]), end = ', ')
        print()


# don't run main if included 
if __name__ == "__main__":
    main()
