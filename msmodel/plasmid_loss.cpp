#include "event_loop.h"

/******************************************************************/
/* DIVERGENCE: SHIFT PLASMID COUNT IN CELL */

std::shared_ptr< V< Cell > > EventLoop::divergeCell_Plasmid(
	E< Cell, Patch, Mlt >& ancestorScope,
	V< Plasmid >& pl,
	int64_t n_events,
	int64_t n_plasmids
) {

	V<Cell>& r = *ancestorScope.source;
	V<Patch>& pch = *ancestorScope.target;
	//auto variant = u.alter_container_shared_ptr(pl, r, n_plasmids);

	/* V2 barcoding */
	Universe::mset< Cell > bc(r);
	Universe::mset_vector< Plasmid > v;
	size_t pl_mult = u.edge< Mlt >(pl, r).data.mult + n_plasmids;
	v.emplace_back(&pl, pl_mult);
	bc.update(v);

	auto variant = u.mset_get_vertex(bc).shared_from_this();

	/* Transition from ancestors to diverged cells. Note that
	n_events, the binomial draw, can be greater than the number of ancestors, hence
	the number of transitions must be capped. */
	int64_t n_transition = std::min((int64_t)ancestorScope.data.mult, n_events);

	if (n_transition != 0) {
		/* Assign diverged cells to patch */
		u.assign(*variant, pch, n_transition);

		/* Remove ancestors from patch */
		u.shift_mult(ancestorScope, -n_transition);


	}

	return variant;
}

/******************************************************************/
/* PLASMID LOSS */

/* Plasmid loss. */
std::shared_ptr< V< Cell > > EventLoop::plasmidLossScoped(
	E<  Plasmid, Cell, Mlt >& plasmidScope,
	E<  Cell, Patch, Mlt >& cellScope
) {
	assert(plasmidScope.target->id == cellScope.source->id);

	V<Plasmid>& pl = *plasmidScope.source;
	V<Cell>& c = *plasmidScope.target;
	V<Patch>& pch = *cellScope.target;

	/* Baseline probability of loss per cell, to be updated based on number of plasmids */
	const double net_loss_prob = pl.data.loss;

	/* Draw n events */
	int64_t n = sampler.rbinom(cellScope.data.mult, net_loss_prob);

	return n == 0 ? nullptr : divergeCell_Plasmid(cellScope, pl, n, -1);
}



/* Plasmid loss for all occurrences of a plasmid */
void EventLoop::plasmidLoss(V< Plasmid >& pl) {

	auto& plScopes = pl.edges_to<Cell, Mlt>();
	for (auto plScope = plScopes.rbegin(); plScope != plScopes.rend(); ++plScope) {
		auto& cellScopes = plScope->target->edges_to<Patch, Mlt>();
		for (auto cellScope = cellScopes.rbegin(); cellScope != cellScopes.rend(); ++cellScope) {
			auto div = plasmidLossScoped(*plScope, *cellScope);
		}
	}

}

/* Global plasmid loss */
void EventLoop::plasmidLossAll() {
	for (auto& pl : u.vertices<Plasmid>()) {
		plasmidLoss(pl);
	}
}