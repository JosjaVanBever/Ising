from exact import get_options
import numpy as np
import numpy.random as random

# shorthand for multiplication of 2 matrices
def matmul(a,b):
    return np.einsum('ij,jk->ik',a,b)

# generate a random symmetric matrix
def rand_sym_mat(D):
    rand_mat = random.rand(D,D)
    return (rand_mat + rand_mat.T) / 2

# output: the orthonormalized eigenvector matrix
def gram_schmidt_columns(X):
    Q, R = np.linalg.qr(X)
    return Q

# average the eigenvalues for general a and b
def enforce_equal_eigenvalues(a,b, symmetric=False):
    # diagonalize a
    eig_a, O_a = np.linalg.eigh(a)  # normalized
    if not symmetric:
        O_a = gram_schmidt_columns(vec_a)
    # diagonalize b
    eig_b, O_b = np.linalg.eigh(b)  # normalized
    if not symmetric:
        O_b = gram_schmidt_columns(vec_b)
    # average the eigenvalues
    eig = (eig_a + eig_b) / 2
    # transform back
    return matmul(np.einsum('ij,j->ij', O_a, eig),O_a.T), \
           matmul(np.einsum('ij,j->ij', O_b, eig),O_b.T)


# class that contains the matrices A(s=up) and
# A(s=down) for a given site s
class site():
    # constructor
    def __init__(self, D):
        # set up and down to real symmetric matrices
        # with equal eigenvalues (spin symmetry)
        self.up, self.down = enforce_equal_eigenvalues(
                rand_sym_mat(D),rand_sym_mat(D), True)
    # get the up (i=0) or down (i=1) matrix
    def get(self,i):
        return [self.up, self.down][i]
    # get a random matrix
    def rand(self):
        return self.get(random.randint(2))
    # representation (for printing)
    def __repr__(self):
        return "Up: %r\nDown:%r\n" \
                % (self.up, self.down)

def main():
    # read the parameters from input
    # N: the amount of sites
    # J, g: parameters of H_{TFIM}
    # D: bond dimension
    N, J, g, D = get_options(True)
    F = 5

    print(N, J, g, D)

    # for each of N sites, keep 2 matrices (for up/down)
    # Remark that these matrices are real symmetric and have common
    # eigenvalues (due to spin symmetry), but does not lead to
    # normalized probabilities p_i: ∑_i p_i != 1.
    # Expectation values should therefore always be normalized.
    site_matrices = [None] + [site(D) for i in range(N)]
    # due to the dummy, site_matrices[1] corresponds to the 1st site

    # site_matrices print as check
    for i in range(1,len(site_matrices)):
        print('Site {0:d}:\n'.format(i), site_matrices[i])

    # the current state of the system is described by an array
    # of bits: 0 = up (= σ^z[0,0]), 1 = down (= σ^z[1,1])
    # initialize at a random state
    state = [None] + [random.randint(2) for i in range(N)]
    # remark that state[1] corresponds to the first site

    # state print as check
    print([state[i] for i in range(1,len(state))])

    # help function to get the current matrix for site i
    def current_site_mat(i):
        return site_matrices[i].get(state[i])

    # RL_matrices[m] = L(m) or R(m) depending on the conditions
    # initialize L(m) = A(s_m), L/R(0) = L/R(N+1) = 1
    RL_matrices = [np.eye(D)] + [current_site_mat(i)
        for i in range(1,len(site_matrices))] + [np.eye(D)]
    # remark that RL_matrices[1] corresponds to site m=1

    # L(m) = A(s_m) * L(m+1), m from N-1 to 1
    for i in range(N,0,-1):  # incl. N, excl. 0
        RL_matrices[i] = matmul(RL_matrices[i], RL_matrices[i+1])

    # remark that W(s) = Tr(A(s_1)...A(s_N)) = Tr(L(1))
    W = np.trace(RL_matrices[1]) # := W(S)

    print('W:', W)

    # estimators for the energy contributions to H_{TFIM}:
    E_z = np.zeros(F)
    E_x = np.zeros(F)
    # estimator for the derivative of ...
    for sweep in range(F):

        # Wp[m] (p='prime') := W(S'_m) = W after flipping spin s_m
        Wp = [None] + [0 for i in range(N)]

        # traverse from s_1 to s_N and do metropolis based flips:
        for m in range(1,len(state)):
            # W(S'_m) = Tr(A(-s_m)*L(m+1)*R(m-1)) GET B(m) HERE!!!
            Aflip = site_matrices[i].get(1-state[i])  # = A(-s_m)
            # Bm = matmul(RL_matrices[m+1], RL_matrices[m-1])
            Wp[m] = np.trace(matmul(Aflip, matmul(RL_matrices[m+1], \
                    RL_matrices[m-1])))
            # change of flipping a spin
            Pflip = min((Wp[m]*Wp[m])/(W*W) ,1)
            # God roles a dice,
            if random.random() < Pflip:
                state[i] = 1 - state[i]  # flips a spin
            # R(m) = R(m-1) * A(s_m)
            RL_matrices[m] = matmul(RL_matrices[m-1], current_site_mat(m))

        print('Wp:', Wp)

        # get the estimator E_z(S) for the diagonal contribution to H_TFIM
        M = 0
        for i in range(1,len(state)-1):
            # true if anti-parallel
            if state[i] ^ state[i+1]:
                M += 1
            else:
                M -= 1
            # periodic boundary conditions
            if state[1] ^ state[N]:
                M += 1
            else:
                M -= 1
        # final result for the estimator E_z
        E_z[sweep] = J * M

        # get the estimator E_x(S) for the non diagonal contribution to H_TFIM
        E_x[sweep] = J * g *sum([Wp[i] for i in range(1,len(Wp))]) / W

        # transverse from s_N to s_1 and evaluate the energy and its derivative:
        for m in range(N,0,-1):
            # L(m) = A(s_m) * L(m+1)
            RL_matrices[m] = matmul(current_site_mat(m), RL_matrices[m+1])



# don't run main if included 
if __name__ == "__main__":
    main()

# RL_matrices = [np.eye(D)] + [mat_i.rand() for mat_i in site_matrices] + [np.eye(D)]