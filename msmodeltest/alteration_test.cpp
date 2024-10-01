#include "pch.h"
#include "utils.h"

#include <string>
#include <algorithm>

TEST(Alteration, alteration) {
	Universe u;
	auto cArch0 = u.insert_vertex_shared_ptr<CellArchetype>(0);
	auto c = u.insertTypedContainerShared(*cArch0, 0);

	auto plArch0 = u.insert_vertex_shared_ptr<PlasmidArchetype>(0);
	for (int i = 4; i >= 0; i--) {
		auto pl_i = u.insertTypedContainerShared(*plArch0, i);
		auto e = u.assign_shared_ptr<Plasmid, Cell>(i, 0);
	}

	/* Alteration: add 1 plasmid */
	/* Plasmid to add */
	auto pl = u.vertex_shared_ptr<Plasmid>(0);
	/* M-shift of +1 */
	int64_t delta = +1;

	/* V2 barcoding: alter barcode before vertex for consistency */
	auto bc = Universe::mset< Cell >(*c);
	bc.update(u.vertex< Plasmid >(0), 1 + delta);
	auto d = u.mset_get_vertex(bc).shared_from_this();
	   
	/* Define an alteration key as an explicit struct */
	Universe::alteration_key< Plasmid, Cell > key(*pl, *c, delta);
	   
	/* Comparison checks */
	ASSERT_EQ(key, (Universe::alteration_key< Plasmid, Cell >(*pl, *c, delta)));
	ASSERT_NE(key, (Universe::alteration_key< Plasmid, Cell >(*pl, *c, delta + 1)));
	ASSERT_NE(key, (Universe::alteration_key< Plasmid, Cell >(*u.vertex_shared_ptr<Plasmid>(1), *c, delta)));

	/* Hash */
	size_t h1 = Universe::alteration_key<Plasmid, Cell>::hash()(Universe::alteration_key<Plasmid, Cell>(*pl, *c, delta));
	size_t h2 = Universe::alteration_key<Plasmid, Cell>::hash()(Universe::alteration_key<Plasmid, Cell>(*pl, *c, delta));
	size_t h3 = Universe::alteration_key<Plasmid, Cell>::hash()(Universe::alteration_key<Plasmid, Cell>(*pl, *c, delta + 1));

	ASSERT_EQ(h1, h2);
	ASSERT_NE(h1, h3);
	   
	/* Register an alteration */
	u.register_alteration(*pl, *c, *d, delta);
	auto& altmap = u.alteration_map<Plasmid, Cell>();

	auto it_f = altmap.find(Universe::alteration_key<Plasmid, Cell>(*pl, *c, delta));
	ASSERT_NE(it_f, altmap.end());
	ASSERT_EQ(it_f->second, &*d);

	auto it_r = altmap.find(Universe::alteration_key<Plasmid, Cell>(*pl, *d, -delta));
	ASSERT_NE(it_r, altmap.end());
	ASSERT_EQ(it_r->second, &*c);

	/* Kill the variant and alterations pointing to and from it */
	u.forget_alteration(*d);

	auto it_f2 = altmap.find(Universe::alteration_key<Plasmid, Cell>(*pl, *c, delta));
	ASSERT_EQ(it_f2, altmap.end());

	auto it_r2 = altmap.find(Universe::alteration_key<Plasmid, Cell>(*pl, *d, -delta));
	ASSERT_EQ(it_r2, altmap.end());
}