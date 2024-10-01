#include "pch.h"

#include "../msmodel/universe.h"

#include "../msgraph/bfs.h"
#include "../msgraph/bfs_counter.h"


/*Replacement NL */

static void setup(Universe& u) {
	u.insert_vertex< Gene >(0);
	
	u.insert_vertex< Chromosome >(0);
	u.assign(u.vertex< Gene >(0), u.vertex< Chromosome >(0), 2);

	u.insert_vertex< Plasmid >(0);
	u.assign(u.vertex< Gene >(0), u.vertex< Plasmid >(0));

	u.insert_vertex< Cell >(0);
	u.assign(u.vertex< Chromosome >(0), u.vertex< Cell >(0));
	u.assign(u.vertex< Plasmid >(0), u.vertex< Cell >(0), 2);

	u.insert_vertex< Patch >(0);
	u.assign(u.vertex< Cell >(0), u.vertex< Patch >(0), 100);
}

TEST(BFS, upwardCount) {
	/* Count instances of an element in all containers */
	Universe u;
	setup(u);

	bfs_callback cb;

	/* Decoupled approach */
	BFS bfs(u);
	bfs.run(u.vertex< Gene >(0), cb);

	/* Count Gene 0 in all containers */
	bfs_counter ct(u);
	bfs.run(u.vertex< Gene >(0), ct);

	/* Check counts */
	EXPECT_EQ(ct.count(u.vertex< Chromosome >(0)), 2);
	EXPECT_EQ(ct.count(u.vertex< Plasmid >(0)), 1);
	EXPECT_EQ(ct.count(u.vertex< Cell >(0)), 4);
	EXPECT_EQ(ct.count(u.vertex< Patch >(0)), 400);
}

TEST(BFS, downwardCount) {
	Universe u;
	setup(u);

	BFS< Universe, bfs_direction::in > bfs(u);
	bfs_counter ct(u);
	bfs.run(u.vertex< Patch >(0), ct);

	/* Check counts */
	EXPECT_EQ(ct.count(u.vertex< Gene >(0)), 400);
	EXPECT_EQ(ct.count(u.vertex< Chromosome >(0)), 100);
	EXPECT_EQ(ct.count(u.vertex< Plasmid >(0)), 200);
	EXPECT_EQ(ct.count(u.vertex< Cell >(0)), 100);
	EXPECT_EQ(ct.count(u.vertex< Patch >(0)), 1);
}
