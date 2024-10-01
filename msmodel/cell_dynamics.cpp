#include "event_loop.h"

/******************************************************************/
/* CELL DIFFUSION */

/* Cell movement along a diffusion edge */
void EventLoop::cellDiffusionScoped(E< Cell, Patch, Mlt >& cellScope, E< Patch, Patch, Diff >& diffScope) {
	/* Consistency check */
	assert(cellScope.target.get() == diffScope.source.get());

	/* Draw no. of moving cells */
	int64_t n = sampler.rbinom(cellScope.data.mult, diffScope.data.diff);

	/* Apply move */
	if (n > 0) {
		u.assign(*cellScope.source, *diffScope.target, +n);
		u.shift_mult(cellScope, -n);
	}
}

/* Cell movement along all relevant diffusion routes */
void EventLoop::cellDiffusion(V< Cell >& c) {

	/* Scope Mlt(Cell, Patch) */
	auto& cellScopes = c.edges_to<Patch, Mlt>();
	for (auto cellScope = cellScopes.rbegin(); cellScope != cellScopes.rend(); ++cellScope) {
		/* Containing patch */
		auto& p = *cellScope->target;

		/* Scope Diff(Patch, Patch) */
		auto& diffScopes = p.edges_to<Patch, Diff>();
		for (auto diffScope = diffScopes.rbegin(); diffScope != diffScopes.rend(); ++diffScope) {
			cellDiffusionScoped(*cellScope, *diffScope);
		}
	}
};

/* Global cell movement */
void EventLoop::cellDiffusionAll() {
	for (auto& it_c : u.vertices<Cell>()) {
		cellDiffusion(it_c);
	}
};

/******************************************************************/
/* CELL BIRTH */

/* Cell reproduction at edge level */
void EventLoop::cellBirthScoped(E<  Cell, Patch, Mlt >& e) {
	Cell& c = e.source->data;
	Patch& p = e.target->data;

	/* Population density in patch */
	const double density = (double)p.popsize / p.capacity;
	/* Abort if full density */
	if (density >= 1.) return;
	/* Damping factor (logistic growth model) */
	const double damp = density > 1 ? 0. : 1. - density;

	/* Binomial draw */
	const double net_reproduction_rate = damp * c.fitness;
	assert(net_reproduction_rate >= 0 && net_reproduction_rate <= 1);
	assert(e.data.mult >= 0 && e.data.mult < 1000000000);
	int64_t n = sampler.rbinom(e.data.mult, net_reproduction_rate);

	if (n != 0) u.shift_mult(e, n);
}


/* Global cell reproduction */
void EventLoop::cellBirthAll() {
	for (auto& it_c : u.vertices<Cell>()) {
		for (auto& cellScope : it_c.edges_to<Patch, Mlt>()) {
			cellBirthScoped(cellScope);
		}
	}
}

/******************************************************************/
/* CELL DEATH */

/* Handle death at edge level */
void EventLoop::cellDeathScoped(E<  Cell, Patch, Mlt >& e) {
	Cell& c = e.source->data;
	Patch& p = e.target->data;

	/* Binomial draw */
	double net_survival = c.survival;
	for (size_t i = 0; i < Cell::n_stress_types; i++) {
		net_survival *= 1. - p.pressure[i] * c.susceptibility[i];
	}
	const double net_death_rate = 1. - net_survival;
	int64_t n = sampler.rbinom(e.data.mult, net_death_rate);

	if (n != 0) u.shift_mult(e, -n);
}

/* Handle global death */
void EventLoop::cellDeathAll() {
	auto& cells = u.vertices<Cell>();
	for (auto cell = cells.rbegin(); cell != cells.rend(); ++cell) {
		auto& cellScopes = cell->edges_to<Patch, Mlt>();
		for (auto cellScope = cellScopes.rbegin(); cellScope != cellScopes.rend(); ++cellScope) {
			cellDeathScoped(*cellScope);
		}
	}
}