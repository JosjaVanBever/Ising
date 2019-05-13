import sys, os 
sys.path.append(os.path.abspath("./include"))
from exact import get_options
from include.mathtools import *
from include.site import Site
import numpy as np
import numpy.random as random


def main():
    # read the parameters from input
    # N: the amount of sites
    # J, g: parameters of H_{TFIM}
    # D: bond dimension
    N, J, g, D = get_options(True)
    F = 5

    print(N, J, g, D)

    # REMARK: site indexing starts at 0 and ends at N-1
    forwards  = range(N)           # 0 ... N-1
    backwards = range(N-1,-1,-1)   # N-1 ... 0

    # for each of N sites, keep 2 matrices (for up/down)
    # Remark that these matrices are real symmetric and have common
    # eigenvalues (due to spin symmetry), but does not lead to
    # normalized probabilities p_i: ∑_i p_i != 1.
    # Expectation values should therefore always be normalized.
    site_matrices = np.array([Site.random(D) for i in forwards])

    # the current state of the system is described by an array
    # of bits: 0 = up (= σ^z[0,0]), 1 = down (= σ^z[1,1])
    # initialize at a random state
    state = [random.randint(2) for i in forwards]

    # help function to get the current matrix for site i
    def current_site_mat(i):
        return site_matrices[i][state[i]]

    # RL_matrices[m] = L(m) or R(m) depending on the conditions
    # Initialize L(m) = A(s_m), L/R(-1) = L/R(N) = 1
    RL_matrices = [current_site_mat(i) for i in forwards] \
                + [np.eye(D)]  # list[-1] == list[len(list)]

    # L(m) = A(s_m) * L(m+1)
    for i in backwards:
        RL_matrices[i] = matmul(RL_matrices[i], RL_matrices[i+1])

    # wave function coefficients for each state S
    W = np.zeros(F)    # W(S) = Tr(A(s_0)...A(s_(N-1)))
    # estimators for the energy contributions to H_{TFIM}:
    E_z = np.zeros(F)  # E_z(S) = J * ∑_{sites} *
                       #    [# anti-parallel(S) - # parallel(S)]
    E_x = np.zeros(F)  # E_x(S) = Jg * ∑_{m} W(S'_m)/W(S)
    # estimator for the derivatives of the energy
    # X_S = ∂W(S)/∂[A(S)]^s_{ij} (= W(S) * [Δ(S)]^s_{ij})
    #     = ∑_{m} [(B(m) + B(m)^t)/2]_{ij} δ_{s_m,s}
    X_S = np.array([Site.zeros(D) for i in forwards])
    # XW = ∑_{S} W(S) * X_S
    XW = np.array([Site.zeros(D) for i in forwards])
    # XWE = ∑_{S} W(S) * E(S) * X_S
    XWE = np.array([Site.zeros(D) for i in forwards])

#    Wderiv = [np.zeros((D,D)) for i in range(F)]   # WRONG DIMENSIONS!!!
    # we should collect matrices per sweep S, per site m, per spin up/down at this site!!!
        # W(S) * Δ(S) = ∂W(S)/∂A(S)
        # <=> W(S) * [Δ(S)]^s_{ij} = ∂W(S)/∂[A(S)]^s_{ij}

    # START EXPERIMENTS
    # site_matrices2 = np.array([site.random(D) for i in forwards])

    print(site_matrices)
    # print(site_matrices2)
    # print(X_S)
    # print('---------------------')

    # site_matrices2 += site_matrices
    # X_S += site_matrices

    # print(site_matrices)
    # print(site_matrices2)
    # print(X_S)
    # print('---------------------')

    site_matrices = 7 * site_matrices

    print(site_matrices)
    # print(site_matrices2)
    # print(X_S)
    # print('---------------------')

    # # X_S[0].set(0,site_matrices[0].get(0))
    # # #X_S[0].up += site_matrices[0].get(0)

    # # print(site_matrices)
    # # print(X_S)

    # #X_S += site_matrices

    # a = np.array([site.zeros(D) for i in forwards])
    # a = a + X_S



    # print(a)

    # X_S[1][0] += site_matrices[1][0]

    # print(a)

    # # site_matrices[1][0] += site_matrices[1][0]

    # # print(site_matrices)
    # # print(X_S)

    # #X_S += site_matrices

    # # print(site_matrices)
    # # print(X_S)

    # eig_a, O_a = np.linalg.eigh(site_matrices[0][0])
    # eig_b, O_b = np.linalg.eigh(site_matrices[0][1])
    # print(eig_a, eig_b)
    # print(matmul(O_a.T,O_a))
    # print(matmul(O_b.T,O_b))



    # do one simulation bin; i.e. collect enough information
    # to calculate the ensemble averages needed to update W(S)
    # each step starts and ends with RL_matrices[m] = L(m)
    for S in range(F):
        # remark that W(s) = Tr(A(s_0)...A(s_(N-1))) = Tr(L(0))
        W[S] = np.trace(RL_matrices[0]) # := W(S)
        # Wp[m] (p='prime') := W(S'_m) = W after flipping spin s_m
        Wp = [0 for i in forwards]
        # traverse from s_0 to s_(N-1) and do metropolis based flips:
        for m in forwards:
            # B(m) = L(m+1)*R(m-1)
            Bm = matmul(RL_matrices[m+1], RL_matrices[m-1])

            # NEW
            X_S[m][state[m]] += (Bm + Bm.T) / 2

            # ∂W(S)/∂[A(S)]^s_{ij} = ∑_{m} [(B(m) + B(m)^t)/2]_{ij}
            # ALARM: DISTINGUISHING UP AND DOWN????????????????????????
            # ASK A(s_m): ar you up or down?
            # when up, then add (Bm + Bm.T) / 2 to (delta_s,s_m in eq. (9)) to Wderiv[sweep S, site m, matrix up]
            # otherwise, add to Wderiv[sweep S, site m, matrix down]
#            Wderiv[S][m][state[m]] += (Bm + Bm.T) / 2

            # get A(-s_m)
            Aflip = site_matrices[m][1-state[m]]
            # W(S'_m) = Tr(A(-s_m)*B(m))
            Wp[m] = np.trace(matmul(Aflip, Bm))
            # change of flipping a spin
            Pflip = min((Wp[m]*Wp[m])/(W[S]*W[S]) ,1)
            # God roles a dice,
            if random.random() < Pflip:
                state[m] = 1 - state[m]  # and flips a spin
            # R(m) = R(m-1) * A(s_m)
            RL_matrices[m] = matmul(RL_matrices[m-1], current_site_mat(m))

        # get the estimator E_z(S) for the diagonal contribution to H_TFIM
        M = 0
        for i in forwards:
            # true if anti-parallel
            if state[i] ^ state[i-1]:  # incl. periodic boundary conditions
                M += 1
            else:
                M -= 1
        # final result for the estimator E_z
        E_z[S] = J * M

        # get the estimator E_x(S) for the non diagonal contribution to H_TFIM
        E_x[S] = J * g * sum(Wp) / W[S]

        # transverse from s_(N-1) to s_0 and evaluate the energy and its derivative:
        for m in backwards:
            # L(m) = A(s_m) * L(m+1)
            RL_matrices[m] = matmul(current_site_mat(m), RL_matrices[m+1])

        XW += W[S] * X_S               # remark that W[S] is a scalar
        XWE += XW * (E_z[S] + E_x[S])  # remark that E_x/z[S] are scalars


    # calculate the ensemble avarages:
    Z = sum(W * W)                       # Z = ∑_{S} W²(S)
    E = (1/Z) * sum(W*W * (E_z + E_x))   # E = <E_z(S) + E_x(S)>_{S}

    # ∂E/∂[A]^s_{ij} = 2 * <[E(S) * [Δ(S)]^s_{ij}]> - 2 * <[Δ(S)]^s_{ij}> * E
    #                = 2/Z * ∑_{S} [W(S) * E(S) * ∂W(S)/∂[A(S)]^s_{ij}]
    #                  - 1/Z * ∑_{S} [W(S) * ∂W(S)/∂[A(S)]^s_{ij}]
    Ederiv = (2/Z) * XWE - (1/Z) * XW

    print('W: ', W)
    print('W*W:', W * W)
    print('E:', E)
    print('Z: ', Z)
    print('Ederiv:', Ederiv)


#     # ∂E/∂[A]^s_{ij} = 2/Z * ∑_{S} [W(S) * E(S) * ∂W(S)/∂[A(S)]^s_{ij}]
#     #                  - 2 * <[Δ(S)]^s_{ij}> * E
# #    Ederiv = (2/Z) * sum(W * (E_z + E_x) * Wderiv) - 2 * Delta * E

#     # [A(S)]^s_{ij} -> [A(S)]^s_{ij} - δ(k) * r^s_{ij} * sgn(∂E/∂[A]^s_{ij})
#     r = rand_sym_mat(D)

#     # <[Δ(S)]^s_{ij}> = 1/Z * ∑_{S} [W(S) * ∂W(S)/∂[A(S)]^s_{ij}]
#     Delta = - 2 (1/Z) * XW




# don't run main if included 
if __name__ == "__main__":
    main()
