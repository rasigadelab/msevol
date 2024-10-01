#include "pch.h"

#include "utils.h"
#include "../msmodel/universe_io.h"

/* Include all relevant V/E types */
static void setup_csv(Universe& u) {

	/* Typed vertices */
	u.insert_vertex< GeneArchetype >(0);
	u.insert_vertex< Gene >(0);
	u.assign< GeneArchetype, Gene >(0, 0);

	u.insert_vertex< ChromosomeArchetype >(0);
	u.insert_vertex< Chromosome >(0);
	u.assign< ChromosomeArchetype, Chromosome >(0, 0);

	u.insert_vertex< PlasmidArchetype >(0);
	u.insert_vertex< Plasmid >(0);
	u.assign< PlasmidArchetype, Plasmid >(0, 0);

	u.insert_vertex< CellArchetype >(0);
	u.insert_vertex< Cell >(0);
	u.assign< CellArchetype, Cell >(0, 0);

	u.insert_vertex< PatchArchetype >(0);
	u.insert_vertex< Patch >(0);
	u.assign< PatchArchetype, Patch >(0, 0);
	u.insert_vertex< Patch >(1);
	u.assign< PatchArchetype, Patch >(0, 1);

	/* Multiplicity edges */
	u.assign< Gene, Plasmid >(0, 0);
	u.assign< Gene, Chromosome >(0, 0);
	u.assign< Plasmid, Cell >(0, 0);
	u.assign< Chromosome, Cell >(0, 0);
	u.assign< Cell, Patch >(0, 0, 100);
	u.assign< Cell, Patch >(0, 1, 200);

	/* Diffusion edges */
	u.insert_edge< Patch, Patch >(0, 1, Diffusion{ 0.01 });
	u.insert_edge< Patch, Patch >(1, 0, Diffusion{ 0.01 });

	/* Plasmid host range edge */
	u.insert_edge< PlasmidArchetype, ChromosomeArchetype, PlasmidRange >(0, 0);
}

TEST(CsvTest, test) {
	constexpr const char* dataDir = "data";

	Universe u;
	setup_csv(u);
	   
	UniverseIO::writeGraph(u, dataDir);
	
	Universe v;
	UniverseIO::readGraph(v, dataDir);
	   
	ASSERT_EQ(u.vertices<GeneArchetype>().size(), v.vertices<GeneArchetype>().size());
	ASSERT_EQ(u.vertices<PlasmidArchetype>().size(), v.vertices<PlasmidArchetype>().size());
	ASSERT_EQ(u.vertices<ChromosomeArchetype>().size(), v.vertices<ChromosomeArchetype>().size());
	ASSERT_EQ(u.vertices<CellArchetype>().size(), v.vertices<CellArchetype>().size());
	ASSERT_EQ(u.vertices<PatchArchetype>().size(), v.vertices<PatchArchetype>().size());

	ASSERT_EQ(u.vertices<Gene>().size(), v.vertices<Gene>().size());
	ASSERT_EQ(u.vertices<Plasmid>().size(), v.vertices<Plasmid>().size());
	ASSERT_EQ(u.vertices<Chromosome>().size(), v.vertices<Chromosome>().size());
	ASSERT_EQ(u.vertices<Cell>().size(), v.vertices<Cell>().size());
	ASSERT_EQ(u.vertices<Patch>().size(), v.vertices<Patch>().size());


	/* Check counts */
	for (size_t i = 0; i < u.vertices<Gene>().size(); i++) {
		ASSERT_EQ(u.vertex_shared_ptr<Gene>(i)->data.count(), v.vertex_shared_ptr<Gene>(i)->data.count());
	}
	for (size_t i = 0; i < u.vertices<Plasmid>().size(); i++) {
		ASSERT_EQ(u.vertex_shared_ptr<Plasmid>(i)->data.count(), v.vertex_shared_ptr<Plasmid>(i)->data.count());
	}
	for (size_t i = 0; i < u.vertices<Chromosome>().size(); i++) {
		ASSERT_EQ(u.vertex_shared_ptr<Chromosome>(i)->data.count(), v.vertex_shared_ptr<Chromosome>(i)->data.count());
	}
	for (size_t i = 0; i < u.vertices<Cell>().size(); i++) {
		ASSERT_EQ(u.vertex_shared_ptr<Cell>(i)->data.count(), v.vertex_shared_ptr<Cell>(i)->data.count());
	}
	for (size_t i = 0; i < u.vertices<Patch>().size(); i++) {
		ASSERT_EQ(u.vertex_shared_ptr<Patch>(i)->data.count(), v.vertex_shared_ptr<Patch>(i)->data.count());
	}

	/* Diffusion edges */
	ASSERT_TRUE((u.has_edge<Patch, Patch, Diff>(0, 1)));
	ASSERT_TRUE((u.has_edge<Patch, Patch, Diff>(1, 0)));

	ASSERT_TRUE((v.has_edge<Patch, Patch, Diff>(0, 1)));
	ASSERT_TRUE((v.has_edge<Patch, Patch, Diff>(1, 0)));

	/* Host range edges */
	ASSERT_TRUE((u.has_edge< PlasmidArchetype, ChromosomeArchetype, PlasmidRange >(0, 0)));
	ASSERT_TRUE((v.has_edge< PlasmidArchetype, ChromosomeArchetype, PlasmidRange >(0, 0)));
}