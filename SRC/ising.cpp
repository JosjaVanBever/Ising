
#include "globals.h"
#include "Htfim.h"
using namespace std;

// get the options from input
int getOptions(int argc, char *argv[], int & L,
			double & J, double & g);

extern"C" {
void dsaupd_(int & IDO, char & BMAT, int & N, char WHICH[2], int & NEV,
		double & TOL, double * RESID, int & NCV, double * V, int & LDV,
		int IPARAM[11], int IPNTR[11], double * WORKD,
		double * WORKL, int & LWORKL, int & INFO );
void dseupd_(bool & RVEC, char & HOWMNY, bool * SELECT, double * D,
			double * Z, int & LDZ, double & SIGMA,
			char & BMAT, int & N, char WHICH[2], int & NEV, double & TOL,
			double * RESID, int & NCV, double * V, int & LDV, int IPARAM[11],
			int IPNTR[11], double * WORKD, double * WORKL, int & LWORKL,
			int & INFO);
}

int main(int argc, char *argv[]) {
	// variables
	int optind;
	int L=-1;
	double J=1;
	double g=0; 

	// get the options
	optind = getOptions(argc, argv, L, J, g);
	// the L option is obliged
	if (L == -1) {
		cerr << syntax << endl;
		abort(); }

	// current state == expansion in basis states!
	// State current(1<<L);

	int ido = 0;			// Reverse communication flag. (in/out-)
	char bmat = 'I';		// Normal eigenvalue equation. (in)
	int n = 1<<L;			// Dimension of the eigenproblem. (in)
	char which[3] = "SM";	// Compute smallest (in magnitude) eigenvalues. (in)
	int nev = 4;			// Number of eigenvalues to be computed. (in)
	double tol = 1e-10;		// Tolerated error on the eigenvalues. (in)

	double resid[n];		// Will contain the final residual vector. (in/out)
	int ncv = min(n, max(2*nev + 1, 20));	// Number of Lanczos vectors. (in)
	double v[n*ncv];		// Will contain the Lanczos vectors. (out)
	int ldv = n;			//???? Leading dimension of v.??? (in)
	int iparam[11];			// Itaration parameters: (in)
		iparam[0] = 1;			// Exact shifts are used.
		iparam[2] = 300;		// Maximum number of Arnoldi iterations allowed.
		iparam[6] = 1;			// Mode 1 of dsaupd (normal eigenvalue equation).
	int ipntr[11];			// Pointers inside working spaces. (work)
	double workd[3*n];		// Workspace that will contain residuals. (out/work)
	int lworkl = ncv*(ncv+8);	// Size of private workspace. (in/work)
	double workl[lworkl];		// Private workspace. (work)
	int info = 0;			// Randomly initial residual vector is used. (in/out)

	bool rvec = true;		// Whether eigenvectors are calculated or not.
	char howmny = 'A';		// Which eigenvectors are calculated. 'A'='All' (in)
	bool select[ncv];		// Eigenvectors that are calculated if howmny = 'S'
	double d[nev];			// Will contain the eigenvalues. (out)
	double sigma;			// Represents the shift if iparam[6] is 3,4,5 (in)
	int ierr = 0;

	assert(nev < ncv);

	dsaupd_(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv,
			iparam, ipntr, workd, workl, lworkl, info );

	// H_TFIM hamiltonian(J,g,L);
	while (ido == 1 | ido == -1) {
		// cout << "1:" << workd[ipntr[0]] << endl;
		// cout << "2:" << workd[ipntr[1]] << endl;
		for (int i=0; i<n; i++) {
			workd[ipntr[1]+i] = 10 * workd[ipntr[0]+i];
		}	
		dsaupd_(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv,
				iparam, ipntr, workd, workl, lworkl, info );
	}

	if (info < 0) {
		cerr << "dsaupd error" << endl;
	} else {
		dseupd_(rvec, howmny, select, d, v, ldv, sigma, 
				bmat, n, which, nev, tol, resid, ncv, v, ldv, 
				iparam, ipntr, workd, workl, lworkl, ierr);
		if (ierr != 0) {
			cout << "dseupd error" << endl;
		} else {
			for (int i=0; i<nev; i++) {
				cout << d[1] << endl;
			}
		}
	}

	// hamiltonian.multiply(workd[ipntr[0]], workd[ipntr[1]]);
	// Bitset chain(L);
	// cout << chain << endl;

	return 0;
}

int getOptions(int argc, char *argv[], int & L,
			double & J, double & g) {
	int c;
	// loop over all options
	while ((c = getopt (argc, argv, "L:J:g:")) != -1) {
	    switch (c) {
		  case 'L':
			L = atoi(optarg); break;
		  case 'J':
			J = atof(optarg); break;
		  case 'g':
			g = atof(optarg); break;
		  case '?':
			cerr << syntax << endl;
			abort();
		  default:
		    abort();
		}
	}
	// return the index of the 1st non-option argument
	return optind;
}

// { "shell_cmd": "cd ..; make clean && make && ./ising.exe -L5 -J1 -g0" }