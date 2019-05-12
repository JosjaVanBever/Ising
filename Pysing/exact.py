# -*- coding: utf-8 -*-
# Program to perform exact diagonalization of the transvers field ising model.

from scipy.sparse.linalg import eigsh
import numpy as np
import getopt, sys

def get_options(getD=False):
    # print the syntax
    def usage():
        program = sys.argv[0].split('/')[-1]
        print("Syntaxis: {0:s} -L <sites> [-J <J> -g <g>]".format(program))
    # default values
    L = -1
    J = 1.0
    g = 0.0
    D = -1
    # get the values from the options
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hL:J:g:D:", ["help"])
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
        elif opt == "-D":
            D = int(arg)
        elif opt in ("-h", "--help"):
            usage()
            sys.exit()
        else:
            assert False, "unhandled option"
    # check obligated parameters
    if L <= 0:
        usage()
        raise IOError("L is obligated and should be larger then 0.")
    if getD and D <= 0:
        usage()
        raise IOError("D is obligated and should be larger then 0.")
    # return the options found
    if getD:
        return L, J, g, D
    else:
        return L, J, g

# H_TFIM = -J * [∑_{<ij>}^L σ^z_i σ^z_j + g * ∑_{i}^L  σ^x_i]
class H_TFIM():
    # constructor
    def __init__(self, L, J, g):
        self.L = L
        self.J = J
        self.g = g
        self.dim = 1 << L
        self.mat = np.zeros((self.dim, self.dim))
        self.mag = np.zeros(self.dim)
        # nearest neighbour term
        for basis in range(self.dim):
            M = 0
            for site in range(L - 1):
                # true if site and site + 1 are anti-parallel
                if (basis >> site ^ basis >> (site + 1)) & 1:
                    M += 1
                else:
                    M -= 1
            # periodic boundary conditions
            # skip this part for open boundary conditions
            if (basis ^ basis >> (L - 1)) & 1:
                M += 1
            else:
                M -= 1
            # factor J * [# anti-parallel - # parallel]
            self.mat[basis, basis] = J * M

        # transvers field term
        for basis in range(self.dim):
            for site in range(L):
                # spinflip at 'site' and factor -Jg
                self.mat[basis ^ (1 << site), basis] = -J * g

        # magnetization
        for basis in range(self.dim):
            for site in range(L):
                # count parallel spins
                if (basis >> site) & 1: 
                    self.mag[basis] += 1
                else:
                    self.mag[basis] -= 1
        # print([self.mag[i] for i in range(len(self.mag))])

    # representation of H_TFIM (for printing)
    def __repr__(self):
        return "-%.1f * [∑_{<ij>}^%d σ^z_i σ^z_j + %.1f * ∑_{i}^%d  σ^x_i]\n%r" \
                % (self.J, self.L, self.g, self.L, self.mat)

def get_order_parameter(vec,ham):
    # vector is supposed to be normalized
    # you can double check by using
    # print(sum([c**2 for c in vec]))
    # the order parameter is defined as:
    #    for a basis vector: |#up - #down| / L
    #    for a general state: weighted in the basis vectors
    return sum([vec[j]**2 * abs(ham.mag[j]) for j in range(len(vec))]) / ham.L

def get_lowest_modes(hamiltonian, k):
    # return the eigenvalues and eigenvectors
    return eigsh(hamiltonian.mat, k=min(k,hamiltonian.dim-1), which='SA')

def analyse_spectrum(L, J, g, print_dominant=False, threshold=0.12):
    # create the hamiltonian
    hamiltonian = H_TFIM(L, J, g)
    # print the hamiltonian
    print('hamiltonian:\n', hamiltonian)

    # diagonalize the hamiltonian
    eigenvalues, eigenvectors = get_lowest_modes(hamiltonian, 30)

    # print the eigenvalues
    print('energies: ', eigenvalues)
    # print the energies per site
    print('ground state energy per site: ', eigenvalues[0] / hamiltonian.L)
    # print the eigenvectors
    print('eigenvector per energy:')
    for i in range(len(eigenvalues)):
        eigenvec = eigenvectors[:, i]
        print('{0:.2f}: '.format(eigenvalues[i]), eigenvec, end=';  ')
        # print out the order parameter for the given state
        print('m = {0:.3f}'.format(get_order_parameter(eigenvec, hamiltonian)))
        # collect the dominant contributions in eigenvec
        if print_dominant:
            dominant = [d for d in range(len(eigenvec)) if eigenvec[d] > threshold]
            for j in range(len(dominant)):
                print('{0:b}'.format(dominant[j]), end=';  ')
            print()

def main():
    # set printing format for numpy
    np.set_printoptions(precision=1, suppress=True, threshold=100)

    # read the parameters from input
    # L: the amount of sites
    # J, g: parameters of H_{TFIM}
    L, J, g = get_options()

    # do a spectrum analysis:
    analyse_spectrum(L,J,g, print_dominant=False)


# don't run main if included 
if __name__ == "__main__":
    main()
