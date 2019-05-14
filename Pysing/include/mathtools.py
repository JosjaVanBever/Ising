# -*- coding: utf-8 -*-
# Extra mathtools used in the ising project.

import numpy as np
import numpy.random as random

# shorthand for multiplication of 2 matrices
# @profile
def matmul(a,b):
    return np.einsum('ij,jk->ik',a,b)

# generate a random symmetric matrix
# with elements in [0,1)
# @profile
def rand_sym_mat(D):
    rand_mat = random.rand(D,D)
    return (rand_mat + rand_mat.T) / 2

# output: the orthonormalized eigenvector matrix
# @profile
def gram_schmidt_columns(X):
    Q, R = np.linalg.qr(X)
    return Q

# average the eigenvalues for general a and b
# @profile
def enforce_equal_eigenvalues(a,b, symmetric=False, normalize=False):
    # diagonalize a
    eig_a, O_a = np.linalg.eigh(a)  # normalized
    if not symmetric:
        O_a = gram_schmidt_columns(O_a)
    # diagonalize b
    eig_b, O_b = np.linalg.eigh(b)  # normalized
    if not symmetric:
        O_b = gram_schmidt_columns(O_b)
    # average the eigenvalues
    eig = (eig_a + eig_b) / 2
    if normalize:
        # normalize the largest eigenvalue to 1
        eig /= eig[0]
    # transform back
    return matmul(np.einsum('ij,j->ij', O_a, eig),O_a.T), \
           matmul(np.einsum('ij,j->ij', O_b, eig),O_b.T)
