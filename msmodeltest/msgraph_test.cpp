#include "pch.h"
#include "utils.h"

TEST(UniverseTest, constructor) {
	Universe u;
}

TEST(UniverseTest, addArchetype) {
	Universe u;

	/* Explicit ID mandatory */
	auto g0 = u.insert_vertex_shared_ptr<GeneArchetype>(0);
	EXPECT_EQ(g0->id, 0);

	/* Pass data */
	auto g1 = u.insert_vertex_shared_ptr(GeneArchetype{ { 0.5, 0.3 }, 0.99 });
	EXPECT_EQ(g1->id, 1);
	EXPECT_EQ(g1->data.susceptibility[0], 0.5);
	EXPECT_EQ(g1->data.susceptibility[1], 0.3);
	EXPECT_EQ(g1->data.fitness, 0.99);
}

TEST(UniverseTest, insertTypedContainerShared) {
	Universe u;
	auto arch = u.insert_vertex_shared_ptr<GeneArchetype>(0);

	/* Explicit ID */
	auto g0 = u.insertTypedContainerShared(*arch, 0);
	EXPECT_EQ(g0->id, 0);
	EXPECT_EQ(g0->data.count(), 0);

	/* Automatic ID */
	auto g1 = u.insertTypedContainerShared(*arch);
	EXPECT_EQ(g1->id, 1);
	EXPECT_EQ(g1->data.count(), 0);

	/* Archetype should have iterable list */
	EXPECT_EQ((arch->edges_to<Gene, Mlt>().size()), u.vertices<Gene>().size());

	/* Archetype should be retrieved from element */
	auto arch_found = u.getArchetypeShared(*g0);
	EXPECT_EQ(arch_found, arch);

	/* Alternative, reference return */
	auto& g2 = u.insertTypedContainer(*arch);
	auto& arch_found2 = u.getArchetype(g2);
	EXPECT_EQ(&arch_found2, &*arch);
}


TEST(UniverseTest, archetypePropagation) {
	Universe u;
	/* Properties of archetype should propagate to element */
	{
		auto a = u.insert_vertex_shared_ptr<GeneArchetype>(0);
		a->data.set(0.1, 0.2);
		for (size_t i = 0; i < pressure_base::n_stress_types; i++) {
			EXPECT_EQ(a->data.susceptibility[i], 0.1);
		}
		EXPECT_EQ(a->data.fitness, 0.2);

		/* Multiple pressure */
		std::array<double, pressure_base::n_stress_types> sv = { 0.3, 0.4 };
		a->data.set(sv, 0.1);
		for (size_t i = 0; i < pressure_base::n_stress_types; i++) {
			EXPECT_EQ(a->data.susceptibility[i], sv[i]);
		}

		auto x = u.insertTypedContainerShared(*a);
		EXPECT_EQ(u.getArchetypeShared(*x), a);
		EXPECT_EQ(x->data.fitness, a->data.fitness);

		for (size_t i = 0; i < pressure_base::n_stress_types; i++) {
			EXPECT_EQ(x->data.susceptibility[i], a->data.susceptibility[i]);
		}
	}
	{
		auto a = u.insert_vertex_shared_ptr<PlasmidArchetype>(0);
		a->data.set(0.1, 0.2, 0.3);
		EXPECT_EQ(a->data.loss, 0.1);
		EXPECT_EQ(a->data.transfer, 0.2);
		EXPECT_EQ(a->data.fitness, 0.3);

		auto x = u.insertTypedContainerShared(*a);
		EXPECT_EQ(u.getArchetypeShared(*x), a);
		EXPECT_EQ(x->data.fitness, a->data.fitness);
	}
	{
		auto a = u.insert_vertex_shared_ptr<ChromosomeArchetype>(0);
		a->data.set(0.1, 0.2);
		EXPECT_EQ(a->data.fitness, 0.1);
		EXPECT_EQ(a->data.survival, 0.2);

		auto x = u.insertTypedContainerShared(*a);
		EXPECT_EQ(u.getArchetypeShared(*x), a);
		EXPECT_EQ(x->data.fitness, a->data.fitness);
		EXPECT_EQ(x->data.survival, a->data.survival);
	}
	{
		auto a = u.insert_vertex_shared_ptr<CellArchetype>(0);
		auto x = u.insertTypedContainerShared(*a);
		EXPECT_EQ(u.getArchetypeShared(*x), a);
	}
	{
		auto a = u.insert_vertex_shared_ptr<PatchArchetype>(0);
		a->data.set(0.1, 0.2);
		EXPECT_EQ(a->data.capacity, 0.1);
		for (size_t i = 0; i < pressure_base::n_stress_types; i++) {
			EXPECT_EQ(a->data.pressure[i], 0.2);
		}

		std::array<double, pressure_base::n_stress_types> sv = { 0.3, 0.4 };
		a->data.set(0.1, sv);
		for (size_t i = 0; i < pressure_base::n_stress_types; i++) {
			EXPECT_EQ(a->data.pressure[i], sv[i]);
		}

		auto x = u.insertTypedContainerShared(*a);
		EXPECT_EQ(u.getArchetypeShared(*x), a);
	}
}

TEST(UniverseTest, destructor) {
	Universe u;
	setup_00(u);
}

TEST(UniverseTest, assign_v2) {
	Universe u;
	setup_00(u);

	/* Explicit multiplicity */
	const size_t m = 100;

	/* Assign should create edge if no edge is present */
	ASSERT_FALSE((u.has_edge< Gene, Chromosome, Mlt >(0, 0)));
	u.assign(u.vertex<Gene>(0), u.vertex<Chromosome>(0), m);
	EXPECT_TRUE((u.has_edge< Gene, Chromosome, Mlt >(0, 0)));
	EXPECT_EQ((u.edge< Gene, Chromosome, Mlt >(0,0).data.mult), m);

	/* Assign should retrieve existing edge if present */
	u.assign(u.vertex<Gene>(0), u.vertex<Chromosome>(0), m);
	EXPECT_EQ((u.edge< Gene, Chromosome, Mlt >(0, 0).data.mult), m * 2);

	/* Assign can find edge by vertex IDs */
	u.assign< Gene, Chromosome >(0, 0, m);
	EXPECT_EQ((u.edge< Gene, Chromosome, Mlt >(0, 0).data.mult), m * 3 );

	/* Assign return edge reference */
	auto& e_ref = u.assign< Gene, Chromosome >(0, 0, m);
	EXPECT_EQ(e_ref.data.mult, m * 4);
}

TEST(UniverseTest, assign_shared_ptr) {
	Universe u;
	setup_00(u);

	/* Explicit multiplicity */
	const size_t m = 100;

	/* Assign should create edge if no edge is present */
	ASSERT_FALSE((u.has_edge< Gene, Chromosome, Mlt >(0, 0)));
	auto e_1 = u.assign_shared_ptr(u.vertex<Gene>(0), u.vertex<Chromosome>(0), m);
	EXPECT_TRUE((u.has_edge< Gene, Chromosome, Mlt >(0, 0)));
	EXPECT_EQ((u.edge< Gene, Chromosome, Mlt >(0, 0).data.mult), m);

	/* Assign should retrieve existing edge if present */
	auto e_2 = u.assign_shared_ptr(u.vertex<Gene>(0), u.vertex<Chromosome>(0), m);
	EXPECT_EQ((u.edge< Gene, Chromosome, Mlt >(0, 0).data.mult), m * 2);

	/* Assign can find edge by vertex IDs */
	auto e_3 = u.assign_shared_ptr< Gene, Chromosome >(0, 0, m);
	EXPECT_EQ((u.edge< Gene, Chromosome, Mlt >(0, 0).data.mult), m * 3);

	/* Assign return edge reference */
	auto e_4 = u.assign_shared_ptr< Gene, Chromosome >(0, 0, m);
	EXPECT_EQ(e_4->data.mult, m * 4);
}

TEST(UniverseTest, shift_count) {
	Universe u;
	setup_00(u);

	/* Gene -x2-> Chrom -x1-> Cell */
	u.assign(u.vertex<Gene>(0), u.vertex<Chromosome>(0), 2);
	u.assign(u.vertex<Gene>(1), u.vertex<Chromosome>(0), 2);
	u.assign(u.vertex<Chromosome>(0), u.vertex<Cell>(0));

	ASSERT_EQ(u.vertex<Cell>(0).data.count(), 0);

	/* Set n(cell) = 10 */
	u.shift_count(u.vertex<Cell>(0), +10);

	EXPECT_EQ(u.vertex<Chromosome>(0).data.count(), 10);
	EXPECT_EQ(u.vertex<Gene>(0).data.count(), 20);
	EXPECT_EQ(u.vertex<Gene>(1).data.count(), 20);
}

TEST(UniverseTest, shift_mult) {
	Universe u;
	setup_00(u);

	/* Gene -x1-> Chromosome */
	u.assign(u.vertex<Gene>(0), u.vertex<Chromosome>(0));
	auto g0 = u.vertex_shared_ptr<Gene>(0);
	auto chr0 = u.vertex_shared_ptr<Chromosome>(0);
	for (size_t i = 0; i < pressure_base::n_stress_types; i++) {
		EXPECT_DOUBLE_EQ(chr0->data.susceptibility[i], g0->data.susceptibility[i]);
	}
	
	EXPECT_DOUBLE_EQ(chr0->data.fitness,
		g0->data.fitness * u.getArchetypeShared(*chr0)->data.fitness);

	for (size_t i = 0; i < pressure_base::n_stress_types; i++) {
		EXPECT_DOUBLE_EQ(chr0->data.susceptibility[i], g0->data.susceptibility[i]);
	}

	/* Add 1 gene */
	{
		auto e = u.edge_shared_ptr<Gene, Chromosome, Mlt>(0, 0);
		u.shift_mult(*e, +1);
		ASSERT_EQ(e->data.mult, 2);
	}
	for (size_t i = 0; i < pressure_base::n_stress_types; i++) {
		EXPECT_DOUBLE_EQ(chr0->data.susceptibility[i], pow(g0->data.susceptibility[i], 2));
	}
	EXPECT_DOUBLE_EQ(chr0->data.fitness, pow(g0->data.fitness, 2) * u.getArchetypeShared(*chr0)->data.fitness);

	for (size_t i = 0; i < pressure_base::n_stress_types; i++) {
		EXPECT_DOUBLE_EQ(chr0->data.susceptibility[i], pow(g0->data.susceptibility[i], 2));
	}
}

TEST(UniverseTest, diffusion) {
	Universe u;
	setup_00(u);

	u.insert_edge(u.vertex<Patch>(0), u.vertex<Patch>(1), Diff{ 0.1 });
	ASSERT_TRUE((u.has_edge<Patch, Patch, Diff>(0, 1)));
}