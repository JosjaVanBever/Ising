
#include "Htfim.h"

void H_TFIM::multiply(const State & in, State & out) {
	// help variables
	int basis, site, M=0;
	// precalculated factors
	double mJg = - J * g, mJ = -J;

	// nearest neighbour term
	for (basis=0; basis<(1<<L); basis++) {
		for (site=0; site<L-1; site++) {
			// check if spins at site and site+1 are parallel
			M += ((basis>>site ^ basis>>(site+1)) & 1)? 1 : -1;
		}
		// periodic boundary conditions
		M += ((basis ^ basis>>L) & 1)? 1 : -1;
		// multiply by -J * [# parallel - # anti-parallel]
		out[basis] = mJ * M * in[basis];
	}
	
	// transverse field term
	for (basis=0; basis<(1<<L); basis++) {
		for (site=0; site<L; site++) {
			// flip the spin at 'site' and multiply by -Jg
			out[basis^(1<<site)] = mJg * in[basis];
		}
	}
}