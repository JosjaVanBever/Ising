import sys, os 
sys.path.append(os.path.abspath("./include"))
from exact import get_options
from include.mathtools import *
from include.site import *
import numpy as np
import numpy.random as random


def main():
    # read the parameters from input
    # N: the amount of sites
    # J, g: parameters of H_{TFIM}
    # D: bond dimension
    N, J, g, D = get_options(True)
    F = F_0 = 80 # e.g. 100
    G = G_0 = 10  # e.g. 10
    delta = delta_0 = 0.05
    Q = 0.95
    k = 1
    tol = 1e-4

    # print out the parameters
    print('INPUT PARAMETERS:\n' + 
        'Number of lattice sites L: %d\n' % (N) +
        'Hamiltonian parameters: J = %.2f, g = %.2f\n' % (J,g) +
        'Bond dimension D = %d\n' % (D) +
        'Initial convergence parameters: ' +
        'δ = %d, G = %d, F = %d\n' % (delta,G,F) +
        'Tolerance on E convergence: %.1e\n' % (tol))
    sys.stdout.flush()  # explicitely flush

    # REMARK: site indexing starts at 0 and ends at N-1
    forwards  = range(N)           # 0 ... N-1
    backwards = range(N-1,-1,-1)   # N-1 ... 0

    # for each of N sites, keep 2 matrices (for up/down)
    # Remark that these matrices are real symmetric and have common
    # eigenvalues (due to spin symmetry), but does not lead to
    # normalized probabilities p_i: ∑_i p_i != 1.
    # Expectation values should therefore always be normalized.
    # These matrices will further be denoted by 'A(S)'.'
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

    # results will be collected in following array
    Energies = []
    StandardDeviations = []

    # is the calculation converged?
    converged = False

    # start convergence:
    print('START CONVERGENCE LOOP:')
    sys.stdout.flush()  # explicitely flush
    while not converged:
        # results are collected per simulation bin
        E_simulation_bin = []

        # Do G simulation bins per set of convergence parameters.
        # Each simulation bin, the ensemble averages are calculated
        # and A(S) is updated.
        for simulation_bin in range(G):

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

            # Do one simulation bin; i.e. collect enough information
            # to calculate the ensemble averages needed to update A(S)
            # After this, do update A(S).
            # Each step of a simulation bin starts and ends
            # with RL_matrices[m] = L(m).
            for S in range(F):
                # remark that W(s) = Tr(A(s_0)...A(s_(N-1))) = Tr(L(0))
                W[S] = np.trace(RL_matrices[0]) # := W(S)
                # Wp[m] (p='prime') := W(S'_m) = W after flipping spin s_m
                Wp = [0 for i in forwards]
                # traverse from s_0 to s_(N-1) and do metropolis based flips:
                for m in forwards:
                    # B(m) = L(m+1)*R(m-1)
                    Bm = matmul(RL_matrices[m+1], RL_matrices[m-1])
                    # ∂W(S)/∂[A(S)]^s_{ij} = ∑_{m} [(B(m) + B(m)^t)/2]_{ij}
                    X_S[m][state[m]] += (Bm + Bm.T) / 2

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

                # transverse from s_(N-1) to s_0 and evaluate the energy
                # and its derivative:
                for m in backwards:
                    # L(m) = A(s_m) * L(m+1)
                    RL_matrices[m] = matmul(current_site_mat(m), RL_matrices[m+1])

                # get estimator E_z(S) for the diagonal contribution to H_TFIM
                M = 0
                for i in forwards:
                    # true if anti-parallel
                    if state[i] ^ state[i-1]:  # incl. periodic boundary conditions
                        M += 1
                    else:
                        M -= 1
                # final result for the estimator E_z
                E_z[S] = J * M

                # get estimator E_x(S) for the non diagonal contribution to H_TFIM
                E_x[S] = J * g * sum(Wp) / W[S]

                # update estimator sums for ∂W(S)/∂[A(S)]^s_{ij} related quantities
                XW += W[S] * X_S               # remark that W[S] is a scalar
                XWE += XW * (E_z[S] + E_x[S])  # remark that E_x/z[S] are scalars

            # calculate the ensemble avarages:
            Z = sum(W * W)                       # Z = ∑_{S} W²(S)
            E = (1/Z) * sum(W*W * (E_z + E_x))   # E = <E_z(S) + E_x(S)>_{S}
            # add the newly evaluated energy to the results
            E_simulation_bin.append(E)

            # ∂E/∂[A]^s_{ij} = 2 * <[E(S)*[Δ(S)]^s_{ij}]> - 2 * <[Δ(S)]^s_{ij}>*E
            #                = 2/Z * ∑_{S} [W(S) * E(S) * ∂W(S)/∂[A(S)]^s_{ij}]
            #                  - 1/Z * ∑_{S} [W(S) * ∂W(S)/∂[A(S)]^s_{ij}]
            Ederiv = (2/Z) * XWE - (1/Z) * XW

            # [A(S)]^s_{ij} -> [A(S)]^s_{ij} - δ(k)*r^s_{ij}*sgn(∂E/∂[A]^s_{ij})
            # COMPLETION OF THIS EQUATION THOROUGHLY VERIFIED, OK!
            r = np.array([Site.random(D) for i in forwards])                          # ALARM COST!!!
            # Actual update of the matrix elements:
            site_matrices -= delta * multiply(r, sign(Ederiv, iterable=True, \
                isSite=True), iterable=True, isSite=True)

            # After each update, spin symmetry is imposed by enforcing equal
            # eigenvalues for the A(s) and A(s) matrix.
            for site in site_matrices:
                site.up, site.down = enforce_equal_eigenvalues(site.up, site.down)

        # after G simulation bins, add the newly collected bin results
        # to ourconvergence results
        Energies.append(np.average(E_simulation_bin))
        StandardDeviations.append(np.std(E_simulation_bin))

        # check wether convergence is reached
        if k > 1:
            if abs(Energies[-1] - Energies[-2]) < tol or k > 3:
                converged = True

        # print current estimations to monitor convergence
        print('k = %d (G = %d, F = %d, δ = %.2f): E = %.8f, σ(E) = %.8f' \
            % (k, G,F,delta, Energies[-1],StandardDeviations[-1]))
        sys.stdout.flush()  # explicitely flush

        # update the convergence parameters and start the G simulation
        # bins all over again
        k += 1       # step monitor
        delta *= Q   # delta = delta_0**k
        F += F_0     # F = k * F_0
        G += G_0     # G = k * G_0
    
    # print out the final result if convergence is reached
    if converged:
        print ('\nConvergence reached!\n' +
            'Final results: E = %.8f, σ(E) = %.8f' \
            % (Energies[-1],StandardDeviations[-1]))


# don't run main if included 
if __name__ == "__main__":
    main()
