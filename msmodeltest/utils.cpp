#include "pch.h"
#include "../msmodel/universe.h"

/* 5 elements of each type, unassigned */
void setup_00(Universe& u) {
	auto& gArch = u.insert_vertex<GeneArchetype>(0);
	gArch.data.set(.5, .8);
	auto& plArch = u.insert_vertex<PlasmidArchetype>(0);
	plArch.data.set(.5, .6, .7);
	auto& chrArch = u.insert_vertex<ChromosomeArchetype>(0);
	chrArch.data.set(.9, .8);
	auto& cllArch = u.insert_vertex<CellArchetype>(0);
	auto& ptchArch = u.insert_vertex<PatchArchetype>(0);
	ptchArch.data.set(1e6, .8);

	for (int i = 0; i < 5; i++) {
		u.insert_vertex< Gene >(i);
		u.assign< GeneArchetype, Gene >(0, i);

		u.insert_vertex< Plasmid >(i);
		u.assign< PlasmidArchetype, Plasmid >(0, i);

		u.insert_vertex< Chromosome >(i);
		u.assign< ChromosomeArchetype, Chromosome >(0, i);

		u.insert_vertex< Cell >(i);
		u.assign< CellArchetype, Cell >(0, i);

		u.insert_vertex< Patch >(i);
		u.assign< PatchArchetype, Patch >(0, i);
	}
}

/* (Gene0 + Gene1) -x2-> Chrom0 -x1-> Cell0 -x100-> Patch0 */
void setup_01(Universe& u) {
	setup_00(u);

	/* Gene -x2-> Chrom -x1-> Cell -x100-> Patch */
	u.assign< Gene, Chromosome >(0, 0, 2);
	u.assign< Gene, Chromosome >(1, 0, 2);
	u.assign< Chromosome, Cell >(0, 0);
	u.assign< Cell, Patch >(0, 0, 100);

	/* Diffusion edges between Patch0, Patch1 */
	u.insert_edge< Patch, Patch, Diff >(0, 1, 0.1);
	u.insert_edge< Patch, Patch, Diff >(1, 0, 0.1);
}

/* Chrom0 -x1-> (Cell x5) -x100-> (Patch x5)*/
void setup_02(Universe& u) {
	auto& chrArch = u.insert_vertex<ChromosomeArchetype>(0);
	chrArch.data.set(.9, .8);
	auto& cllArch = u.insert_vertex<CellArchetype>(0);
	auto& ptchArch = u.insert_vertex<PatchArchetype>(0);
	ptchArch.data.set(1e6, .8);

	for (int i = 0; i < 5; i++) {
		u.insert_vertex< Chromosome >(i);
		u.assign< ChromosomeArchetype, Chromosome >(0, i);

		u.insert_vertex< Cell >(i);
		u.assign< CellArchetype, Cell >(0, i);

		u.insert_vertex< Patch >(i);
		u.assign< PatchArchetype, Patch >(0, i);
	}

	for (int i = 0; i < 5; i++) {
		u.assign<Cell, Patch>(i, i, 100);
	}
}

/* Plasmid */