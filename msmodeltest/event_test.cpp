#include "pch.h"
#include "utils.h"
#include "../msmodel/event_loop.h"
using namespace std;



TEST(EventLoop, constructor) {
	Universe u;
	setup_01(u);
	Sampler sp;
	EventLoop ev(u, sp);
}



TEST(EventLoop, cellBirthScoped) {
	Universe u;
	setup_01(u);
	Sampler sp;
	EventLoop ev(u, sp);

	auto e = u.edge_shared_ptr<Cell, Patch, Mlt>(0, 0);
	size_t multBefore = e->data.mult;

	ev.cellBirthScoped(*e);
	ASSERT_GT(e->data.mult, multBefore);
}

TEST(EventLoop, cellBirthAll) {
	Universe u;
	setup_02(u); // 5 cells, each x100 in one of 5 patches
	Sampler sp;
	EventLoop ev(u, sp);

	vector<size_t> multBefore;
	for (size_t i = 0; i < u.vertices<Cell>().size(); i++) {
		auto e = u.edge_shared_ptr<Cell, Patch, Mlt>(i, i);
		multBefore.push_back(e->data.mult);
	}

	ev.cellBirthAll();

	for (size_t i = 0; i < u.vertices<Cell>().size(); i++) {
		auto e = u.edge_shared_ptr<Cell, Patch, Mlt>(i, i);
		ASSERT_GT(e->data.mult, multBefore[i]);
	}
}



TEST(EventLoop, cellDeathScoped) {
	Universe u;
	setup_01(u);
	Sampler sp;
	EventLoop ev(u, sp);

	auto e = u.edge_shared_ptr<Cell, Patch, Mlt>(0, 0);
	size_t multBefore = e->data.mult;

	ev.cellDeathScoped(*e);
	ASSERT_LT(e->data.mult, multBefore);
}



TEST(EventLoop, cellDeathAll) {
	Universe u;
	setup_02(u); // 5 cells, each x100 in one of 5 patches
	Sampler sp;
	EventLoop ev(u, sp);

	vector<size_t> multBefore;
	for (size_t i = 0; i < u.vertices<Cell>().size(); i++) {
		auto e = u.edge_shared_ptr<Cell, Patch, Mlt>(i, i);
		multBefore.push_back(e->data.mult);
	}

	ev.cellDeathAll();

	for (size_t i = 0; i < u.vertices<Cell>().size(); i++) {
		auto e = u.edge_shared_ptr<Cell, Patch, Mlt>(i, i);
		ASSERT_LT(e->data.mult, multBefore[i]);
	}
}



TEST(EventLoop, plasmidLossScoped) {
	Universe u;
	setup_00(u);
	Sampler sp;
	EventLoop ev(u, sp);

	/* Assign plasmid */
	{ auto e = u.assign_shared_ptr<Gene, Plasmid>(0, 0); }
	{ auto e = u.assign_shared_ptr<Plasmid, Cell>(0, 0); }
	{ auto e = u.assign_shared_ptr<Cell, Patch>(0, 0, 100); }

	/* Scope edges */
	auto plasmidScope = u.edge_shared_ptr<Plasmid, Cell, Mlt>(0, 0);
	auto cellScope = u.edge_shared_ptr<Cell, Patch, Mlt>(0, 0);

	/* Deferred divergence pointer */
	auto div = ev.plasmidLossScoped(*plasmidScope, *cellScope);

	ASSERT_FALSE(div == nullptr);
	ASSERT_TRUE(div->data.count() > 0);

	Cell& c = cellScope->source->data;
	ASSERT_TRUE(div->data.fitness > c.fitness);
	ASSERT_TRUE(div->data.survival == c.survival);
	ASSERT_DOUBLE_EQ(div->data.fitness, c.fitness / plasmidScope->source->data.fitness);
}



TEST(EventLoop, plasmidLoss) {
	Universe u;
	setup_00(u);
	u.mset_register_all();
	Sampler sp;
	EventLoop ev(u, sp);

	/* Assign plasmid in several cells */
	auto e = u.assign_shared_ptr<Gene, Plasmid>(0, 0);
	size_t nCells_before = u.vertices<Cell>().size();
	for (size_t i = 0; i < nCells_before; i++) {
		auto ep = u.assign_shared_ptr<Plasmid, Cell>(0, i);
		auto ec = u.assign_shared_ptr<Cell, Patch>(i, 0, 100);
	}

	/* Plasmid loss in all cells */
	auto pl = u.vertex_shared_ptr<Plasmid>(0);
	ev.plasmidLoss(*pl);

	size_t nCells_after = u.vertices<Cell>().size();

	/* CMS tracking: no new cell is created because cells without plasmid preexist */
	ASSERT_EQ(nCells_after, nCells_before);
}

TEST(EventLoop, plasmidLossAll) {
	Universe u;
	setup_00(u);
	u.mset_register_all();
	Sampler sp;
	EventLoop ev(u, sp);

	size_t nPlasmids = u.vertices<Plasmid>().size();
	size_t nCells = u.vertices<Cell>().size();
	for (size_t i = 0; i < nCells; i++) {
		for (size_t j = 0; j < nPlasmids; j++) {
			auto ep = u.assign_shared_ptr<Plasmid, Cell>(j, i);
		}
		auto ec = u.assign_shared_ptr<Cell, Patch>(i, 0, 100);
	}

	ev.plasmidLossAll();

	size_t nCells_after = u.vertices<Cell>().size();
	ASSERT_GT(nCells_after, nCells * nPlasmids);
}

TEST(EventLoop, cellDiffusionScoped) {
	Universe u;
	setup_00(u);
	u.mset_register_all();
	Sampler sp;
	EventLoop ev(u, sp);

	/* Set diffusion, unidirectional */
	u.insert_edge<Patch, Patch, Diff>(0, 1, 0.1);

	/* Assign cells to patch 0 */
	auto e = u.assign_shared_ptr<Cell, Patch>(0, 0, 1000);

	/* Scope: Mlt(Cell, Patch) and Diff(Patch, Patch) */
	auto cellScope = u.edge_shared_ptr<Cell, Patch, Mlt>(0, 0);
	auto diffScope = u.edge_shared_ptr<Patch, Patch, Diff>(0, 1);

	ev.cellDiffusionScoped(*cellScope, *diffScope);

	ASSERT_LT((u.edge_shared_ptr<Cell, Patch, Mlt>(0, 0)->data.mult), 1000);
	ASSERT_GT((u.edge_shared_ptr<Cell, Patch, Mlt>(0, 1)->data.mult), 0);
}

TEST(EventLoop, cellDiffusion) {
	Universe u;
	setup_00(u);
	u.mset_register_all();
	Sampler sp;
	EventLoop ev(u, sp);

	/* Set diffusion, unidirectional */
	u.insert_edge<Patch, Patch, Diff>(0, 1, 0.1);
	u.insert_edge<Patch, Patch, Diff>(2, 3, 0.1);

	/* Assign cells to patches 0 and 2 */
	{ auto e = u.assign_shared_ptr<Cell, Patch>(0, 0, 1000); }
	{ auto e = u.assign_shared_ptr<Cell, Patch>(0, 2, 1000); }

	/* Cell should diffuse for patch0 to patch1 and from patch2 to patch3 */
	ev.cellDiffusion(*u.vertex_shared_ptr<Cell>(0));

	ASSERT_LT((u.edge_shared_ptr<Cell, Patch, Mlt>(0, 0)->data.mult), 1000);
	ASSERT_GT((u.edge_shared_ptr<Cell, Patch, Mlt>(0, 1)->data.mult), 0);

	ASSERT_LT((u.edge_shared_ptr<Cell, Patch, Mlt>(0, 2)->data.mult), 1000);
	ASSERT_GT((u.edge_shared_ptr<Cell, Patch, Mlt>(0, 3)->data.mult), 0);
}

TEST(EventLoop, cellDiffusionAll) {
	Universe u;
	setup_00(u);
	Sampler sp;
	EventLoop ev(u, sp);

	/* Set diffusion, unidirectional */
	u.insert_edge<Patch, Patch, Diff>(0, 1, 0.1);
	u.insert_edge<Patch, Patch, Diff>(2, 3, 0.1);

	/* Assign cells to patches 0 and 2 */
	{ auto e = u.assign_shared_ptr<Cell, Patch>(0, 0, 1000); }
	{ auto e = u.assign_shared_ptr<Cell, Patch>(0, 2, 1000); }
	{ auto e = u.assign_shared_ptr<Cell, Patch>(1, 0, 1000); }
	{ auto e = u.assign_shared_ptr<Cell, Patch>(1, 2, 1000); }

	/* Cell should diffuse for patch0 to patch1 and from patch2 to patch3 */
	ev.cellDiffusionAll();

	ASSERT_LT((u.edge_shared_ptr<Cell, Patch, Mlt>(0, 0)->data.mult), 1000);
	ASSERT_GT((u.edge_shared_ptr<Cell, Patch, Mlt>(0, 1)->data.mult), 0);
	ASSERT_LT((u.edge_shared_ptr<Cell, Patch, Mlt>(0, 2)->data.mult), 1000);
	ASSERT_GT((u.edge_shared_ptr<Cell, Patch, Mlt>(0, 3)->data.mult), 0);

	ASSERT_LT((u.edge_shared_ptr<Cell, Patch, Mlt>(1, 0)->data.mult), 1000);
	ASSERT_GT((u.edge_shared_ptr<Cell, Patch, Mlt>(1, 1)->data.mult), 0);
	ASSERT_LT((u.edge_shared_ptr<Cell, Patch, Mlt>(1, 2)->data.mult), 1000);
	ASSERT_GT((u.edge_shared_ptr<Cell, Patch, Mlt>(1, 3)->data.mult), 0);
}