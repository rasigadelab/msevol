#pragma once
#include <random>
#include <vector>
#include <unordered_map>
#include "universe.h"
#include "../msgraph/sampler.h"

/* Clarify vertex typename */
template<typename _VData>
using V = Universe::template V< _VData >;

/* Clarify edge typename */
template<typename _Source, typename _Target, typename _EData>
using E = Universe::template E< _Source, _Target, _EData >;

class EventLoop {
protected:
	Sampler& sampler;
public:
	Universe& u;
	EventLoop(Universe& u_in, Sampler& sampler) : u{ u_in }, sampler{ sampler } {};

	/* Random number generator */
	std::default_random_engine rng;

	/* Single binomial draw */
	inline int64_t rbinom(const size_t n, const double prob) {
		std::binomial_distribution<int64_t> rbinom_(n, prob);
		return rbinom_(rng);
	}
	   
protected:
	std::shared_ptr< V< Cell > > divergeCell_Plasmid(
		E< Cell, Patch, Mlt >& ancestorScope,
		V< Plasmid >& pl,
		int64_t n_events,
		int64_t n_plasmids
	);

public:
	/* Cell movement along a diffusion edge */
	void cellDiffusionScoped(E<Cell, Patch, Mlt>& cellScope, E<Patch, Patch, Diff>& diffScope);

	/* Cell movement along all relevant diffusion routes */
	void cellDiffusion(V< Cell >& c);

	/* Global cell movement */
	void cellDiffusionAll();

	/* Cell reproduction at edge level */
	void cellBirthScoped(E< Cell, Patch, Mlt >& e);
	/* Global cell reproduction */
	void cellBirthAll();

	/* Cell death at edge level */
	void cellDeathScoped(E< Cell, Patch, Mlt >& e);
	/* Global cell death */
	void cellDeathAll();

/**************************************************************************************/
/* PLASMID LOSS */

	/* Plasmid loss. Parameter div is the "deferred divergence" pointer,
	can be null. Specify scopes using edge parameters. Will fail if run
	twice using same scope and diverged cell.
	*/
	std::shared_ptr< V< Cell > > plasmidLossScoped(
		E< Plasmid, Cell, Mlt >& plasmidScope,
		E< Cell, Patch, Mlt >& cellScope
	);

	/* Scope overload: plasmid loss for all occurrences of a plasmid */
	void plasmidLoss(V< Plasmid >& pl);

	/* Global plasmid loss */
	void plasmidLossAll();

/**************************************************************************************/
/* PLASMID TRANSFER */
	/* in plasmid_transfer.cpp */
public:
	/* Count no. of plasmids with a given archetype in a cell */
	size_t plasmidCountInCell(V< PlasmidArchetype >& pl_a, V< Cell >& cell);
public:
	/* Compute no. of donors and plasmids in each patch */
	void updateConjugationEdges(V< Plasmid >& plasmid);

	/* Compute all plasmid transfer events */
	void plasmidTransferAll();

	/**************************************************************************************/
	/* RUN ALL EVENTS */
public:
	/* */

	/* Run all events once */
	void step();
	void randomStep();

	/* Run all events several times */
	void run(const size_t duration);
};