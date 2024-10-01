/****************************************************************************************************/
/* PLASMID TRANSFER EVENTS */
/*

REQUIRED STRUCTURES:
CONFIGURATION:
	- plasmid host range is materialized by edges with Range type between a PlasmidArchetype and a
	ChromosomeArchetype. Plasmids can transfer to any cell whose chromosome is in range.
	- PlasmidMovement::transfer is the nominal probability of transfer when the no. of donors equals
	patch capacity. The net transfer probability is the nominal probability
	times the density of recipients in patch.
	- PlasmidMovement::maxcount is the maximum no. of plasmids of the same archetype that a cell can
	harbor. Use maxcount == 1 to simulate plasmid exclusion. Larger values can be used to control
	model complexity. Plasmids that attempt to transfer to a saturated cell are lost.

BOOKKEEPING:
	- the model keeps track of the no. of transferable plasmids and donors
	in a patch in Edge< Plasmid, Patch, Conjugation >

STRATEGY:

The implementation avoids checking transfer between all cell types, which has complexity
quadratic in the no. of cell populations in a patch (which can grow quickly). To avoid
this O(n²) time, we purposedly ignore which cell transfers the plasmid to which recipient.

To prevent a single cell to gain several identical plasmids per step, the recipient
cells that receive a plasmid (= the variants) are removed from the pool of recipients 
by moving them to a temporary, helper patch. The variants are not added immediately to the patch.
Their number is stored in a separate structure. The variants are added to their
respective patches once all transfer events have been performed for a given plasmid.

The same cell can receive several plasmids of different types during the same step.

DETAILS:
- for each plasmid type, compute no. of potential recipients in each patch. These are the cells whose
	chromosome is in host range of the plasmid type and that harbor less plasmids of this
	type than the maximum allowed. Store result in E< PlasmidArchetype, Patch, Recipients >.
- for each plasmid in each patch:
	- compute no. of plasmids, store results in	E< Plasmid, Patch, Conjugation >
	- compute no. of donors for the mass-action model
	- distribute the transferred plasmids across recipients
*/
/****************************************************************************************************/
#include <iostream>
#include "event_loop.h"
using std::cout;
using std::endl;

/* Show detailed debug info ? */
//#define DEBUG_INFO

/* Helper: count no. of plasmids with a given archetype in a cell */
size_t EventLoop::plasmidCountInCell(V< PlasmidArchetype >& pl_a, V< Cell >& cell) {
	size_t plasmid_count = 0;
	for (auto& e_pl : cell.edges_from< Plasmid, Mlt >()) {
		if (u.has_edge< Mlt >(pl_a, *e_pl.source)) {
			plasmid_count += e_pl.data.mult;
		}
	}
	return plasmid_count;
}

void EventLoop::updateConjugationEdges(V< Plasmid >& plasmid) {

	/* Plasmid types */
	V< PlasmidArchetype >& pl_a = u.getArchetype(plasmid);

	/* Reset counters of existing edges */
	for (auto& e_conj : plasmid.edges_to< Patch, Conjugation >()) {
		e_conj.data.n_donors = 0;
		e_conj.data.n_plasmids = 0;
	}

	/* Compute no. of plasmids and donors in each patch */
	for (auto& e_cell : plasmid.edges_to< Cell, Mlt >()) {

		auto& cell = *e_cell.target;

		/* Scan patches */
		for (auto& e_patch : cell.edges_to< Patch, Mlt >()) {

			auto& patch = *e_patch.target;

			/* Retrieve of create conjugation edge */
			auto e_conjugation = u.edge_shared_ptr< Conjugation >(plasmid, patch);
			if(!e_conjugation) e_conjugation = u.insert_edge_shared_ptr< Conjugation >(plasmid, patch);

			/* Update donor and plasmid count in patch */
			assert(e_conjugation != nullptr);
			e_conjugation->data.n_donors += e_patch.data.mult;
			e_conjugation->data.n_plasmids += e_cell.data.mult * e_patch.data.mult;
		}
	}
}

void EventLoop::plasmidTransferAll() {

	/* Scan plasmid types */
	for (auto& pl_a : u.vertices< PlasmidArchetype >()) {

		/* Scan plasmids */
		for (auto& e_pl : pl_a.edges_to< Plasmid, Mlt >()) {
			auto& plasmid = *e_pl.target;

			if (plasmid.data.count() == 0) continue;

			/* Compute no. of waiting plasmids in each patch */
			updateConjugationEdges(plasmid);

			/* 2-step approach: store no. of variants to be added in each patch after plasmid uptake,
			then update multiplicity after all cells and patches have been scanned. This is to avoid
			multiple plasmid acquisition by the same cell during the step. */
			struct transition {
				std::shared_ptr< V< Cell > > cell;
				std::shared_ptr< V< Patch > > patch;
				int64_t n;
				transition(
					const std::shared_ptr < V< Cell > >& cell,
					const std::shared_ptr < V< Patch > >& patch,
					const int64_t n
				) : cell{ cell }, patch{ patch }, n{ n } {}
			};
			std::vector< transition > transitions;

			/* Temporary patch to retain altered cells */
			auto& tmp_container = u.insert_vertex< Patch >();

			/* Scan compatible chromosome archetypes */
			for (auto& e_chrArch : pl_a.edges_to< ChromosomeArchetype, PlasmidRange >()) {

				/* Scan compatible chromosomes */
				for (auto& e_chr : e_chrArch.target->edges_to< Chromosome, Mlt >()) {

					/* Scan compatible cells */
					for (auto& e_cell : e_chr.target->edges_to< Cell, Mlt >()) {

						auto& cell = *e_cell.target;

						/* Check whether cell can accept more plasmids of the same type */
						if (plasmidCountInCell(pl_a, cell) >= pl_a.data.max_count) continue;

						/* Deferred divergence pointer: perform alteration iif the event occurs */
						std::shared_ptr< V< Cell > > variant_ptr = nullptr;

						/* Scan patches */
						for (auto& e_patch : cell.edges_to< Patch, Mlt >()) {
							auto& patch = *e_patch.target;

							if (&patch == &tmp_container) continue;

							/* Conjugation edge holds bookkeeping info */
							auto e_conjugation = u.edge_shared_ptr< Conjugation >(plasmid, patch);

							/* FIXME temporary hack
							 * If patch contains potential recipient but no donor the conjugation edge is not created.
                              Abort loop in this case */

							if (!e_conjugation || (e_conjugation->data.n_plasmids == 0)) continue;

							/* Net uptake rate, mass-action model:	(1 - (1 - rho)^m) d / K
							where rho = transfer rate, m = no. of plasmids per donor, d / K = donor density
							*/
							double m = double(e_conjugation->data.n_plasmids) / double(e_conjugation->data.n_donors);

							double rate = (1. - std::pow(1. - plasmid.data.transfer, m)) *
								double(e_conjugation->data.n_donors) / double(patch.data.capacity);

							if (rate > 1.) rate = 1.;

							/* No. of transfers */
							int64_t n = sampler.rbinom(e_patch.data.mult, rate);

							/* DEBUG INFO */
#ifdef DEBUG_INFO
							cout << "Plasmid #" << plasmid.id
								<< ", Cell #" << cell.id
								<< ", Patch #" << patch.id
								<< ", mult = " << e_patch.data.mult
								<< ", n_plasmids = " << e_conjugation->data.n_plasmids
								<< ", rate = " << rate
								<< ", n_events = " << n
								<< endl;
#endif

							if (n > 0) {
								/* Realize alteration if needed */
								if (variant_ptr == nullptr) {

									/* V2 barcoding */
									Universe::mset< Cell > bc(cell);
									Universe::mset_vector< Plasmid > v;

									size_t pl_mult = 1;
									auto e_plasmid_cell = u.edge_ptr< Mlt >(plasmid, cell);
									if (e_plasmid_cell) pl_mult += e_plasmid_cell->data.mult;

									v.emplace_back(&plasmid, pl_mult);
									bc.update(v);

									variant_ptr = u.mset_get_vertex(bc).shared_from_this();
								}

								/* Move recipients to the temporary container to avoid
								their disappearance. */
								u.assign(cell, tmp_container, n);
								u.shift_mult(e_patch, -n);

								/* Register transition. Cells that have received a plasmid
								are temporarily excluded from the recipients until the end of the
								uptake step. This is to avoid that the same cell receives several
								plasmids during step. */
								transitions.emplace_back(std::move(variant_ptr), std::move(patch.shared_from_this()), n);
							} // END TRANSITION
						} // PATCH SCAN LOOP
					} // END COMPATIBLE CELLS
				} // END COMPATIBLE CHROMOSOMES
			} // END COMPATIBLE CHROMOSOME ARCHETYPES

			/* Step 2: apply transitions and change multiplicities. */
			for (auto& it : transitions) {
				u.assign(*it.cell, *it.patch, it.n);
			}

			/* Remove temporary recipients */
			for (auto& it : tmp_container.edges_from< Cell, Mlt >()) {
				u.shift_mult(it, -int64_t(it.data.mult));
			}
			u.erase_vertex(tmp_container);

		} // END PLASMID
	} // END plasmidArchetype
}