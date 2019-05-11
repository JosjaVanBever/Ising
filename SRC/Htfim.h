
#ifndef HTFIM
#define HTFIM

#include "globals.h"
typedef std::vector<double> State;

// H_TFIM = -J * [∑_{ij} σ^z_i σ^z_j + g * ∑_{i}  σ^x_i]
class H_TFIM {
	public:
		// constructor
		H_TFIM(double J, double g, int L) : J(J), g(g), L(L) {};
		// out = H_TFIM * in
		void multiply(const State & in, State & out);
	private:
		// parameters in H_TFIM
		double J, g;
		// parameters of the lattice
		int L;
};

#endif