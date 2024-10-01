#include "pch.h"

#include "pch.h"
#include "utils.h"
#include "../msmodel/event_loop.h"
using namespace std;

using std::cout;
using std::endl;

TEST(EventLoop, plasmidTransfer) {
	Universe u;
	Sampler sp;
	EventLoop ev(u, sp);

	u.insert_vertex(0, PlasmidArchetype(0., 0.01, 1000, 1.));
	u.insertTypedContainer(u.vertex< PlasmidArchetype >(0), 0);
	u.insertTypedContainer(u.vertex< PlasmidArchetype >(0), 1);

	u.insert_vertex(0, ChromosomeArchetype(0., 1.));
	u.insertTypedContainer(u.vertex< ChromosomeArchetype >(0), 0);

	u.insert_vertex(0, CellArchetype());
	u.insertTypedContainer(u.vertex< CellArchetype >(0), 0);

	u.insert_vertex(0, PatchArchetype(1000, {0., 0.}));
	u.insertTypedContainer(u.vertex< PatchArchetype >(0), 0);
	u.insertTypedContainer(u.vertex< PatchArchetype >(0), 1);

	u.assign(u.vertex< Plasmid >(0), u.vertex< Cell >(0));
	u.assign(u.vertex< Plasmid >(1), u.vertex< Cell >(0));
	u.assign(u.vertex< Chromosome >(0), u.vertex< Cell >(0));
	u.assign(u.vertex< Cell >(0), u.vertex< Patch >(0), 1000);
	u.assign(u.vertex< Cell >(0), u.vertex< Patch >(1),  500);

	u.insert_edge(u.vertex< PlasmidArchetype >(0), u.vertex< ChromosomeArchetype >(0), PlasmidRange());

#ifdef VERBOSE
	cout << endl;
	for (size_t i = 0; i < 10; ++i) {
		cout << "Iter #" << i << endl;
		ev.plasmidTransferAll();
	}
#endif
}