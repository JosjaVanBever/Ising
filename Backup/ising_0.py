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
class H_TFIM(np.ndarray):
    # constructor
    def __init__(self, L, J, g):
        self._L = L
        self._J = J
        self._g = g
        self.dtype = np.dtype(np.float16)
        self.shape = (L,L) # to be changed!

    # representation of H_TFIM (for printing)
    def __repr__(self):
        return "-%.1f * [∑_{<ij>}^%d σ^z_i σ^z_j + %.1f * ∑_{i}^%d  σ^x_i]" \
                % (self._J, self._L, self._g, self._L)

    def multiply(self, input):
        return input
        # do nothing

def main():    
    # read the parameters from input
    # L: the amount of sites
    # J, g: parameters of H_{TFIM}
    L, J, g = get_options()

    print(L, J, g)
    hamiltonian = H_TFIM(L, J, g)
    print(hamiltonian)

    identity = np.eye(8)
    print(identity)
    test = np.array([0, 5, 6, 8, 7, 9, 4, 5])
    print(identity * test[0])
    #print(identity * test)
    #print(np.dot(identity,test))
    #print(np.einsum('ij,j->i', identity, test))
    eigenvalues, eigenvectors = eigsh(hamiltonian, k=6)
    #print(eigenvalues)
    #print(eigenvectors.shape)

if __name__ == "__main__":
    main()

# customized build
# { "shell_cmd": "python3 ising.py -L8 -J1 -g0.5"}