#include "pch.h"

#include "../hgraph/hgraph.h"
#include "minigraph.h"

TEST(hgraph, constructor) {
	Graph g;
}

TEST(hgraph, insert_vertex) {
	Graph g;

	/* Explicit ID */
	auto a0 = g.insert_vertex_shared_ptr<A>(0);
	ASSERT_EQ(a0->id, 0);
	ASSERT_TRUE(a0->active());

	/* Automatic ID */
	auto a1 = g.insert_vertex_shared_ptr<A>();
	ASSERT_EQ(a1->id, 1);

	auto b0 = g.insert_vertex_shared_ptr<B>();
	ASSERT_EQ(b0->id, 0);

	/* Vertex maps increment */
	ASSERT_EQ(g.vertices<A>().size(), 2);
	ASSERT_EQ(g.vertices<B>().size(), 1);

	/* Check iteration */
	size_t A_count = 0;
	for(auto& vertex : g.vertices<A>()) {
		A_count++;
	}
	ASSERT_EQ(g.vertices<A>().size(), A_count);

	size_t B_count = 0;
	for (auto& vertex : g.vertices<B>()) {
		B_count++;
	}
	ASSERT_EQ(g.vertices<B>().size(), B_count);

	/* Additional return types */
	Graph::V<A>& a_ref = g.insert_vertex<A>();
	Graph::V<A>* a_ptr = g.insert_vertex_ptr<A>();
	std::shared_ptr< Graph::V<A> > a_shared = g.insert_vertex_shared_ptr<A>();
}

//#ifdef WIP__




TEST(hgraph, vertexLookup) {
	Graph g;
	auto a0 = g.insert_vertex_shared_ptr<A>(0);
	auto b0 = g.insert_vertex_shared_ptr<B>(0);

	ASSERT_EQ(g.vertex_shared_ptr<A>(0), a0);
	ASSERT_EQ(g.vertex_shared_ptr<B>(0), b0);

	ASSERT_EQ(g.vertex_shared_ptr<A>(1), nullptr);
	ASSERT_EQ(g.vertex_shared_ptr<B>(1), nullptr);

	ASSERT_TRUE(g.has_vertex<A>(0));
	ASSERT_TRUE(g.has_vertex<B>(0));

	ASSERT_FALSE(g.has_vertex<A>(1));
	ASSERT_FALSE(g.has_vertex<B>(1));

	/* Additional return types */
	Graph::V<A>& a_ref = g.vertex<A>(0);
	Graph::V<A>* a_ptr = g.vertex_ptr<A>(0);
	std::shared_ptr< Graph::V<A> > a_shared = g.vertex_shared_ptr<A>(0);

	Graph::V<A>* a_ptr_orig = a0.get();

	ASSERT_EQ(a_shared.get(), a_ptr_orig);
	ASSERT_EQ(&a_ref, a_ptr_orig);
	ASSERT_EQ(a_ptr, a_ptr_orig);

	/* vertex should throw */
	ASSERT_THROW(g.vertex<A>(999), std::out_of_range);

	/* Check static_assert fails */
	//g.vertex< double >(0);
}

TEST(hgraph, vertexLookup2) {
	Graph g;
	auto& a = g.insert_vertex< A >(0);

	Graph::V<A>& a_ref = g.vertex<A>(0);
	Graph::V<A>* a_ptr = g.vertex_ptr<A>(0);
	std::shared_ptr< Graph::V<A> > a_shared = g.vertex_shared_ptr<A>(0);

	EXPECT_EQ(&a_ref, &a);
	EXPECT_EQ(a_ptr, &a);
	EXPECT_EQ(a_shared.get(), &a);
	EXPECT_EQ(&*a_shared, &a);
}

TEST(hgraph, insertEdgeByReference) {
	Graph g;
	auto a0 = g.insert_vertex_shared_ptr<A>(0);
	auto b0 = g.insert_vertex_shared_ptr<B>(0);

	auto e_a0b0 = g.insert_edge_shared_ptr(*a0, *b0, X{});
	ASSERT_TRUE(e_a0b0->active());
	g.insert_edge(*a0, *b0);
	//g.insert_edge_shared_ptr(b0, a0); // Undeclared

	ASSERT_NE((g.vertex_shared_ptr<A>(0)->edges_to<B, X>().front_ptr()), nullptr);
	ASSERT_NE((g.vertex_shared_ptr<B>(0)->edges_from<A, X>().front_ptr()), nullptr);

	/* Additional return types */
	Graph::E< A, B >& e_ref = g.insert_edge( g.insert_vertex<A>(), g.insert_vertex<B>());
	Graph::E< A, B >* e_ptr = g.insert_edge_ptr(g.insert_vertex<A>(), g.insert_vertex<B>());

	/* Additional syntax for edge tags without data */
	g.insert_edge< X >( g.insert_vertex<A>(), g.insert_vertex<B>() );
}

TEST(hgraph, insertEdgeByID) {
	Graph g;
	auto a0 = g.insert_vertex_shared_ptr<A>(0);
	auto b0 = g.insert_vertex_shared_ptr<B>(0);

	g.insert_edge<A, B, X>(0, 0);
	g.insert_edge<A, B>(0, 0);

	ASSERT_NE((g.vertex_shared_ptr<A>(0)->edges_to<B, X>().front_ptr()), nullptr);
	ASSERT_NE((g.vertex_shared_ptr<B>(0)->edges_from<A, X>().front_ptr()), nullptr);
}

TEST(hgraph, edgeLookup) {
	Graph g;
	g.insert_vertex<A>(0);
	g.insert_vertex<A>(1);
	g.insert_vertex<B>(0);
	g.insert_vertex<B>(1);

	auto ref1 = g.insert_edge_shared_ptr<A, B, X>(0, 0);
	auto ref2 = g.insert_edge_shared_ptr<A, B>(0, 0);

	auto source = g.vertex_shared_ptr<A>(0);
	auto target = g.vertex_shared_ptr<B>(0);
	auto badtarget = g.vertex_shared_ptr<B>(1);

	/* Find by ID (returns pointer) */
	ASSERT_EQ(ref1, (g.edge_shared_ptr<A, B, X>(0, 0)));
	ASSERT_EQ(nullptr, (g.edge_shared_ptr<A, B, X>(0, 1)));
	ASSERT_EQ(ref2, (g.edge_shared_ptr<A, B>(0, 0)));
	ASSERT_EQ(nullptr, (g.edge_shared_ptr<A, B>(0, 1)));

	/* Find by reference (returns pointer) */
	ASSERT_EQ(ref1, (g.edge_shared_ptr<X>(*source, *target)));
	ASSERT_EQ(nullptr, (g.edge_shared_ptr<X>(*source, *badtarget)));
	ASSERT_EQ(ref2, (g.edge_shared_ptr<>(*source, *target)));
	ASSERT_EQ(nullptr, (g.edge_shared_ptr<>(*source, *badtarget)));

	/* Has by ID (returns bool) */
	ASSERT_TRUE((g.has_edge<A, B, X>(0, 0)));
	ASSERT_FALSE((g.has_edge<A, B, X>(0, 1)));
	ASSERT_TRUE((g.has_edge<A, B>(0, 0)));
	ASSERT_FALSE((g.has_edge<A, B>(0, 1)));

	/* Has by reference (returns bool) */
	ASSERT_TRUE((g.has_edge<X>(*source, *target)));
	ASSERT_FALSE((g.has_edge<X>(*source, *badtarget)));
	ASSERT_TRUE((g.has_edge<>(*source, *target)));
	ASSERT_FALSE((g.has_edge<>(*source, *badtarget)));

	/* Additional return types */
	Graph::E<A, B>& e_ref = g.edge(*source, *target);
	Graph::E<A, B>* e_ptr = g.edge_ptr(*source, *target);
}

TEST(hgraph, erase_edge) {
	Graph g;
	g.insert_vertex<A>(0);
	g.insert_vertex<B>(0);
	g.insert_edge(*g.vertex_shared_ptr<A>(0), *g.vertex_shared_ptr<B>(0), X{});

	g.erase_edge(*g.edge_shared_ptr<A, B, X>(0, 0));
	ASSERT_FALSE((g.has_edge<A, B, X>(0, 0)));

	ASSERT_EQ((g.vertex_shared_ptr<A>(0)->edges_to<B, X>().front_ptr()), nullptr);
	ASSERT_EQ((g.vertex_shared_ptr<B>(0)->edges_from<A, X>().front_ptr()), nullptr);
}

TEST(hgraph, erase_vertex) {
	Graph g;
	g.insert_vertex<A>(0);
	g.insert_vertex<B>(0);
	g.insert_edge(*g.vertex_shared_ptr<A>(0), *g.vertex_shared_ptr<B>(0), X{});

	g.erase_vertex(*g.vertex_shared_ptr<A>(0));

	/* Check vertex removal in graph map and list */
	ASSERT_FALSE(g.has_vertex<A>(0));
	ASSERT_EQ((g.vertices<A>().front_ptr()), nullptr);

	/* Check edge removal */
	ASSERT_FALSE((g.has_edge<A, B, X>(0, 0)));
	ASSERT_EQ((g.vertex_shared_ptr<B>(0)->edges_from<A, X>().front_ptr()), nullptr);
}

TEST(hgraph, foreach_edge_type) {
	Graph g;
	g.insert_vertex<A>(0);
	g.insert_vertex<B>(0);
	g.insert_vertex<B>(1);
	g.insert_vertex<B>(2);
	g.insert_edge(g.vertex<A>(0), g.vertex<B>(0), X{});
	g.insert_edge(g.vertex<A>(0), g.vertex<B>(1), Y{});
	g.insert_edge(g.vertex<A>(0), g.vertex<B>(2), XY{});

	/* Selecting X should traverse A0->B0 and A0->B2 */
	size_t counter{ 0 };
	g.foreach_edge_type<X>(g.vertex_shared_ptr<A>(0)->edges_out_tuple, [&](auto& x) { counter++; });
	EXPECT_EQ(counter, 2);

	/* Selecting Y should traverse A0->B1 and A0->B2 */
	counter = 0;
	g.foreach_edge_type<Y>(g.vertex_shared_ptr<A>(0)->edges_out_tuple, [&](auto& x) { counter++; });
	EXPECT_EQ(counter, 2);

	/* Selecting XY should traverse only A0->B2 */
	counter = 0;
	g.foreach_edge_type<XY>(g.vertex_shared_ptr<A>(0)->edges_out_tuple, [&](auto& x) { counter++; });
	EXPECT_EQ(counter, 1);
}

TEST(hgraph, counters) {
	/* Count all vertices and edges and return a pair */
	Graph g;
	g.insert_vertex<A>(0);
	g.insert_vertex<B>(0);
	g.insert_vertex<B>(1);
	g.insert_vertex<B>(2);
	g.insert_edge(g.vertex<A>(0), g.vertex<B>(0), X{});
	g.insert_edge(g.vertex<A>(0), g.vertex<B>(1), Y{});
	g.insert_edge(g.vertex<A>(0), g.vertex<B>(2), XY{});

	size_t count = g.order();
	ASSERT_EQ(count, 4);

	auto size = g.size();
	ASSERT_EQ(size.first, 4);
	ASSERT_EQ(size.second, 3);

	/* Specific counters */
}

TEST(hgraph, iterateVertices) {
	Graph g;

	const size_t n = 5;
	for (size_t i = 0; i < n; ++i) {
		g.insert_vertex<A>(i);
	}

	auto it_begin = g.vertices< A >().begin();
	ASSERT_EQ((*it_begin).id, 0); // Operator*
	ASSERT_EQ(it_begin->id, 0); // Operator->
	ASSERT_FALSE(it_begin.is_null());

	++it_begin;
	ASSERT_EQ((*it_begin).id, 1); // Operator*
	ASSERT_EQ(it_begin->id, 1); // Operator->
	ASSERT_FALSE(it_begin.is_null());

	auto it_end = g.vertices< A >().end();
	ASSERT_TRUE(it_end.is_null());

	auto& list = g.vertices<A>();

	size_t i_test = 0;
	for (auto& it : list) {
		ASSERT_EQ(it.id, i_test);
		i_test++;
	}
}

TEST(hgraph, noncopyable) {
	/* Code below SHOULD NOT COMPILE */
	Graph g;

//#define SHOULD_NOT_COMPILE

#ifdef SHOULD_NOT_COMPILE
	/* Avoid nasty bug due to copying instead of reference */
	auto a = g.insert_vertex< A >();
	auto b = g.insert_vertex< B >();
	auto e = g.insert_edge< A, B >(0, 0);
#endif

	for (size_t i = 0; i < 5; ++i) {
		g.insert_vertex< A >(i);
	}

	/* Reference-based range-for is OK */

	size_t i_test = 0;
	for (auto& it : g.vertices< A >()) {
		ASSERT_EQ(it.id, i_test);
		i_test++;
	}

#ifdef SHOULD_NOT_COMPILE
	/* Copy-based range-for is NOT OK */
	i_test = 0;
	for (auto it : g.vertices< A >()) {
		ASSERT_EQ(it.id, i_test);
		i_test++;
	}
#endif

}

TEST(hgraph, vertexLifecycle) {
	/* Delete vertices during iteration and check safety */
	Graph g;

	g.insert_vertex< A >();

	/* Check refcounts and iterators */
	auto a0 = g.vertex_shared_ptr< A >(0);
	ASSERT_EQ(a0.use_count(), 2);

	/* Remove from graph but keep alive */
	g.erase_vertex(*a0);
	ASSERT_EQ(a0.use_count(), 1);
	ASSERT_FALSE(a0->active());
	ASSERT_EQ(g.vertices< A >().size(), 1); // List should still contain "dying" vertex

	/* Release shared_ptr */
	a0.reset();
	ASSERT_EQ(g.vertices< A >().size(), 0); // List should not contain vertex


	/* Try again with iterator */
	g.insert_vertex< A >();
	ASSERT_EQ(g.vertices< A >().size(), 1);

	auto it0 = g.vertices< A >().begin();
	ASSERT_EQ(it0.use_count(), 2);

	g.erase_vertex(*it0);
	ASSERT_EQ(g.vertices< A >().size(), 1); // List should still contain "dying" vertex

	/* Release iterator */
	it0.reset();
	ASSERT_EQ(g.vertices< A >().size(), 0); // List should not contain vertex
}

TEST(hgraph, vertexIterateLifecycle) {
	Graph g;

	const size_t n = 5;
	for (size_t i = 0; i < n; ++i) {
		g.insert_vertex<A>(i);
	}

	size_t i_test = 0;
	for(auto it = g.vertices< A >().begin(); it != g.vertices< A >().end(); ++it) {
		/* Call deletion */
		g.erase_vertex(*it);
		/* Check that deletion was delayed */
		ASSERT_EQ(it->id, i_test);
		/* Check that vertex is marked as inactive */
		ASSERT_FALSE(it->active());
		i_test++;
	}

	/* Check that all vertices are properly removed after iteration */
	ASSERT_EQ(g.vertices< A >().size(), 0);
}

TEST(hgraph, vertexRangeFor) {
	/* Vertices are deleted twice when using range for but not standard for loop */
	/* Only difference it iterator dereference when */

	Graph g;

	const size_t n = 5;
	for (size_t i = 0; i < n; ++i) {
		g.insert_vertex<A>(i);
	}

	size_t i_test = 0;
	for (auto& it : g.vertices< A >()) {
		/* Call deletion */
		g.erase_vertex(it);
		/* Check that deletion was delayed */
		ASSERT_EQ(it.id, i_test);
		/* Check that vertex is marked as inactive */
		ASSERT_FALSE(it.active());
		i_test++;
	}

	/* Check that all vertices are properly removed after iteration */
	ASSERT_EQ(g.vertices< A >().size(), 0);
}

TEST(hgraph, edgeLifecycle) {
	Graph g;

	g.insert_vertex<A>(0);
	g.insert_vertex<B>(0);
	g.insert_edge(*g.vertex_shared_ptr<A>(0), *g.vertex_shared_ptr<B>(0));

	/* Locked edge should prevent vertex deletion */
	auto e = g.edge_shared_ptr<A, B>(0, 0);
	ASSERT_EQ(e.use_count(), 2);

	g.erase_vertex(*g.vertex_shared_ptr<A>(0));
	ASSERT_EQ(e.use_count(), 1);
	ASSERT_EQ(g.vertices< A >().size(), 1);

	/* Edge should not be active anymore */
	ASSERT_FALSE((g.has_edge< A, B >(0, 0)));

	/* Unlocking edge should allow vertex destruction */
	e.reset();
	ASSERT_EQ(g.vertices< A >().size(), 0);


}

struct ShuffleCfg {
	using vertex_types = std::tuple< A, B >;
	using edge_types = std::tuple <
		Edge<A, B>, Edge< B, A >
	>;
};

class ShuffleGraph : public graph< ShuffleCfg > {};

TEST(hgraph, shuffle) {

	/* Shuffle order of vertices and edges in all lists */

	ShuffleGraph g;

	size_t n = 3;

	for (size_t i = 0; i < n; i++) {
		g.insert_vertex< A >(i);
		g.insert_vertex< B >(i);
		g.insert_edge(*g.vertex_shared_ptr< A >(i), *g.vertex_shared_ptr< B >(i));
		g.insert_edge(*g.vertex_shared_ptr< B >(i), *g.vertex_shared_ptr< A >(i));
	}

	g.shuffle();

	/* FIXME should use ASSERT_EQ on order */

	//for (auto& it : g.vertices< A >()) {
	//	std::cout << it.id << " ";
	//}
	//std::cout << std::endl;
}
//#endif // WIP__