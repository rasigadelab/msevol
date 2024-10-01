#pragma once

/****************************************************************************/
/* EDGE PROPERTIES */

/* As per msgraph_io requirements, all properties *must* expose a constructor
compatible with the make_from_tuple construction defined in universe_io.h */

/* Diffusion constant governs (passive) movement of containees
between containers */
struct Diffusion {
	double diff;
	Diffusion(const double diff = 0.) : diff{ diff } {}
};
using Diff = Diffusion;

/* Host range compatibility between plasmid and chromosome archetypes */
struct PlasmidRange {};

/* Plasmid transfer bookkeeping */
struct Conjugation {
	/* Number of donor cells in patch */
	size_t n_donors;
	/* Total number of plasmids in patch */
	size_t n_plasmids;
	Conjugation() : n_donors{ 0 }, n_plasmids{ 0 } {}
};

