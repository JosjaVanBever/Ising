# -*- coding: utf-8 -*-
# Program to perform calculations for the transvers field ising model.

from scipy.sparse.linalg import eigsh
import numpy as np
import random
import sys
import getopt, sys

def usage():
    "Syntaxis: ising.exe -L <sites> [-J <J> -g <g>]"

def get_options():
    # default values
    L = -1
    J = 1.0
    g = 0.0
    # get the values from the options
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hL:J:g:", ["help"])
    except getopt.GetoptError as err:
        # print help information and exit:
        print(str(err))  # will print sth like "option -a not recognized"
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-L":
            L = int(arg)
        elif opt == "-J":
            J = float(arg)
        elif opt == "-g":
            g = float(arg)
        elif opt in ("-h", "--help"):
            usage()
            sys.exit()
        else:
            assert False, "unhandled option"
    # check obligated parameters
    if L <= 0:
        raise IOError("L is obligated and should be larger then 0.")
        usage()
    # return the options found
    return L, J, g

# H_TFIM = -J * [∑_{<ij>}^L σ^z_i σ^z_j + g * ∑_{i}^L  σ^x_i]
class H_TFIM():
    # constructor
    def __init__(self, L, J, g):
        self.L = L
        self.J = J
        self.g = g
        size = 1 << L
        self.mat = np.zeros((size, size))
        # nearest neighbour term
        for basis in range(size):
            M = 0
            for site in range(L - 1):
                # true if site and site + 1 are anti-parallel
                if (basis >> site ^ basis >> (site + 1)) & 1:
                    M += 1
                else:
                    M -= 1
            # periodic boundary conditions
            if (basis ^ basis >> L) & 1:
                M += 1
            else:
                M -= 1
            # factor J * [# anti-parallel - # parallel]
            self.mat[basis, basis] = J * M

        # transvers field term
        for basis in range(size):
            for site in range(L - 1):
                # spinflip at 'site' and factor -Jg
                self.mat[basis ^ (1 << site), basis] = -J * g

    # representation of H_TFIM (for printing)
    def __repr__(self):
        return "-%.1f * [∑_{<ij>}^%d σ^z_i σ^z_j + %.1f * ∑_{i}^%d  σ^x_i]\n%r" \
                % (self.J, self.L, self.g, self.L, self.mat)

def main():    
    # read the parameters from input
    # L: the amount of sites
    # J, g: parameters of H_{TFIM}
    L, J, g = get_options()

    print(L, J, g)
    hamiltonian = H_TFIM(L, J, g)
    print(hamiltonian)

    size = 1 << L
    ham = np.zeros((size, size))
    # nearest neighbour term
    for basis in range(size):
        M = 0
        for site in range(L - 1):
            # true if site and site + 1 are anti-parallel
            if (basis >> site ^ basis >> (site + 1)) & 1:
                M += 1
            else:
                M -= 1
        # periodic boundary conditions
        if (basis ^ basis >> L) & 1:
            M += 1
        else:
            M -= 1
        # factor J * [# anti-parallel - # parallel]
        ham[basis, basis] = J * M

    # transvers field term
    for basis in range(size):
        for site in range(L - 1):
            # spinflip at 'site' and factor -Jg
            ham[basis ^ (1 << site), basis] = -J * g

    print(ham)

    identity = np.eye(8)
    #print(identity)
    #test = np.array([0, 5, 6, 8, 7, 9, 4, 5])
    #print(identity * test[0])
    #print(identity * test)
    #print(np.dot(identity,test))
    #print(np.einsum('ij,j->i', identity, test))
    eigenvalues, eigenvectors = eigsh(hamiltonian.mat, k=6)
    print(eigenvalues)
    #print(eigenvectors.shape)

if __name__ == "__main__":
    main()

# customized build
# { "shell_cmd": "python3 ising.py -L8 -J1 -g0.5"}


# print("check: " + str(3 >> 1))
#     print("check: " + str(3 >> 0))
#     print("check: " + str(1 ^ 3))
#     print("check: " + str(2 & 1))
#     print("check: " + str(0 >> 0 ^ 0 >> 1))
#     print("check: " + str((0 >> 0 ^ 0 >> 1) & 1))