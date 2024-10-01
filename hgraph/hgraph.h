#pragma once

#include <tuple>
#include <array>	
#include <unordered_map>
#include <memory> // std::shared_ptr
#include <variant> // std::monostate
#include <stdexcept>
#include <cassert>
#include <random>
#include <algorithm>

#include "variadic.h"
#include "mlist.h"

/*****************************************************************************************/
/* HETEROGENEOUS TYPE-SAFE MULTIGRAPH */

/*
MOTIVATION.
	In-memory graph container with compile-time definition of allowed vertex and edge types.
	Built for high-frequency lookup/insertion/removal rather than iteration.
	 
FEATURES.
	- non-intrusive container: vertices and edges optionally hold arbitrary data types in a common ::data member.
		Data types need not be modified to be attached to a graph.
	- strong schema enforcement: allowed types of vertices and edges define a schema that is enforced
		during compilation. Attempts to insert unexpected vertices or edges will not compile.
	- constant-time lookup/rewriting: vertex/edge retrieval, insertion and deletion is amortized 
		O(1) independent of graph size or vertex degree
	- safe iteration: iterators over vertices and edges is never invalidated even when another thread or
		function removes the current vertex or edge from the graph.
	- compile-time traversal rules: subgraph selection based on vertex/edge labels uses metaprogramming
		with zero runtime overhead.

BASIC USAGE.
	- pass a schema configuration struct as template parameters, containing aliases vertex_types and edge_types
	- declare vertex types as a tuple of typenames and edge types as a tuple of edge signatures, eg,
		struct my_schema {
			using vertex_types = std::tuple< A, B, C >;
			using edge_types = std::tuple< Edge< A, B >, Edge< A, C > >;
		}
	- instantiate by passing vertex and edge types to the template constructor,
		hgraph::graph< my_schema > myGraph

EDGE LIFECYCLE.
	- each edge E is managed by a shared_ptr owned by the graph.
	- default edge iterators or lookup return shared_ptr's
	- calling graph::erase_edge removes E from the graph, but deletion is deferred
		until all other shared_ptr's are released
	- between removal from graph and effective destruction, E is in 'inactive but reachable' state:
		E.active() returns false; graph.has_edge(E) returns false; however E is still
		reachable by any iterator created prior to E's removal

VERTEX LIFECYCLE.
	- each vertex V is managed by a shared_ptr owned by the graph and optional shared_ptr's
		owned by edges pointing to V
	- calling graph::erase_vertex removes V and all its associated edges from the graph; deletion
		is deferred until all shared_ptr's are released, including iterators or edges pointing to V
	- between removal from graph and final release, V is in 'inactive but reachable' state:
		V.active() returns false; graph.has_vertex(V) returns false; however V is still
		reachable by any iterator created prior to V's removal

*/

namespace hgraph {

	/* Edge placeholder */
	template< typename _Source, typename _Target, typename _EData = std::monostate>
	struct Edge {
		/* Public access to types */
		using source_t = _Source;
		using target_t = _Target;
		using data_t = _EData;
	};

	/* Hashmap key for edges */
	using edge_map_key = std::tuple<size_t, size_t>;

/***********************************************************************/
/*
Type-safe directed multigraph.
Allowed vertex types should be passed as a tuple of types.
Allowed edge types should be passed as a tuple of Edge< > triplets.
*/
	template< typename _Cfg >
	class graph
	{
	public:
		/* Access to declared vertex types */
		using vertex_types = typename _Cfg::vertex_types;
		/* Access to declared edge types */
		using edge_types = typename _Cfg::edge_types;

		/* Vertex forward declaration */
		template<typename _VData = std::monostate> class V;
		/* Edge forward declaration */
		template<typename _Source, typename _Target, typename _EData = std::monostate> class E;

/***********************************************************************/
/* VALIDITY CHECK */
	private:
		/* Static check for vertex validity */
		template<typename _VData>
		struct assert_vertex_is_valid {
			static_assert(
				vrdc::has_type< _VData, vertex_types >::value,
				"Undeclared vertex data type."
				);
		};


		/* Static check for edge validity */
		template<typename _Source, typename _Target, typename _EData>
		struct edge_is_valid_base {
			static_assert(
				vrdc::has_type< _Source, vertex_types >::value,
				"Undeclared source data type in edge."
				);

			static_assert(
				vrdc::has_type< _Target, vertex_types >::value,
				"Undeclared target data type in edge."
				);

			static_assert(
				vrdc::has_type<Edge< _Source, _Target, _EData >, edge_types >::value,
				"Undeclared edge type.");
		};

		/* Static check for edge validity */
		template<typename _Source, typename _Target, typename _EData>
		struct assert_edge_is_valid : edge_is_valid_base< _Source, _Target, _EData > {};

		/* Static check for edge validity */
		template<typename _Edge>
		using assert_edgetype_is_valid = edge_is_valid_base<
			typename _Edge::source_t,
			typename _Edge::target_t,
			typename _Edge::data_t
		>;

/***********************************************************************/
/* EDGE DEFINITION */
	public:

		/* Type-safe edge.*/
		template<typename _Source, typename _Target, typename _EData>
		class E :
			public mlist_hook< E< _Source, _Target, _EData >, 2 >,
			public Edge< _Source, _Target, _EData >,
			assert_edge_is_valid< _Source, _Target, _EData >
		{
			template<typename> friend class graph;
		public:

			/* Source vertex of the edge */
			std::shared_ptr< V< _Source > > source;
			/* Target vertex of the edge */
			std::shared_ptr< V< _Target > > target;
			/* Optional data field */
			_EData data;
			/* Edge is active ?*/
			inline bool active() { return active_; };

			~E() {
				/* Remove edge from source's outgoing edges and target's incoming edges */
				source->edges_to< _Target, _EData >().erase(*this);
				target->edges_from< _Source, _EData >().erase(*this);
			}

			/* Edge is non-copyable */
			E(const E&) = delete;
			E& operator=(const E&) = delete;

		protected:
			/* Protected constructor */
			explicit E(
				V< _Source >& _source,
				V< _Target >& _target,
				_EData _data = {}) :
				source{ _source.shared_from_this() }, target{ _target.shared_from_this() }, data{ _data }, active_{ true } {}

			/* Protected constructor */
			explicit E(
				std::shared_ptr< V< _Source > > _source,
				std::shared_ptr< V< _Target > > _target,
				_EData _data = {}) :
				source{ _source }, target{ _target }, data{ _data }, active_{ true } {}

			/* Flag: edge is active in graph */
			bool active_;
		};

		/* Transform an edge declaration into a proper hgraph edge */
		template<typename _Edge>
		using as_hgraph_edge = E< typename _Edge::source_t, typename _Edge::target_t, typename _Edge::data_t >;


/***********************************************************************/
/* VERTEX DEFINITION */
	public:

		/* Type-safe vertex. */
		template<typename _VData>
		class V :
			public mlist_hook< V<_VData> >,
			assert_vertex_is_valid<_VData>
		{
			template<typename> friend class graph;
		public:

			/* Unique vertex identifier */
			const size_t id;

			/* Optional vertex data field */
			using data_t = _VData;
			_VData data;

			/* edge is active ?*/
			inline bool active() { return active_; };

			//private:
				/* Tuple of outgoing edges, multilist handle = 0 by convention */
			template<typename... T> struct as_edges_out {
				using tuple_t = std::tuple< as_hgraph_edge< T >... >;
				using ilist_t = std::tuple< mlist< as_hgraph_edge< T >, 0>... >;
			};

			/* Tuple of incoming edges, multilist handle = 1 by convention */
			template<typename... T> struct as_edges_in {
				using tuple_t = std::tuple< as_hgraph_edge< T >... >;
				using ilist_t = std::tuple< mlist< as_hgraph_edge< T >, 1>... >;
			};

			/* Predictate for selection of edges */
			template<typename _Edge>
			struct edge_has_source : std::is_same<typename _Edge::source_t, _VData> {};

			/* Predictate for selection of edges */
			template<typename _Edge>
			struct edge_has_target : std::is_same<typename _Edge::target_t, _VData> {};

			/* Filter valid edge types for this vertex */
			using edges_out_tuple_t = typename vrdc::repack_if< edge_has_source, edge_types >::type;
			using edges_in_tuple_t = typename vrdc::repack_if< edge_has_target, edge_types >::type;

			/* Build edge lists */
			using edges_out_t = typename vrdc::repack_as< as_edges_out, edges_out_tuple_t >::type;
			using edges_in_t = typename vrdc::repack_as< as_edges_in, edges_in_tuple_t >::type;

		public:
			/* Install edge lists */
			typename edges_out_t::ilist_t edges_out_tuple;
			typename edges_in_t::ilist_t edges_in_tuple;

		public:
			/* Type-based access to first and last edge */
			template< typename _Source, typename _EData = std::monostate>
			inline auto& edges_from() {
				using edge_t = E<_Source, _VData, _EData>;
				return std::get< mlist<edge_t, 1> >(edges_in_tuple);
			}

			/* Type-based access to the map of outgoing edges with specific child and data types */
			template< typename _Target, typename _EData = std::monostate>
			inline auto& edges_to() {
				using edge_t = E<_VData, _Target, _EData>;
				return std::get< mlist<edge_t, 0> >(edges_out_tuple);
			}

			~V() {
				graph_.vertices<_VData>().erase(*this);
			}

			/* Vertex is non-copyable */
			V(const V&) = delete;
			V& operator=(const V&) = delete;

		protected:
			/* Protected constructor */
			explicit V(graph< _Cfg >& _graph, const size_t _id, _VData _data = {}) :
				graph_{ _graph }, id{ _id }, data{ _data }, active_{ true } {}
			/* Flag: vertex is active in graph */
			bool active_;
			/* Reference to parent graph, required to call vertices in destructor */
			graph< _Cfg >& graph_;
		};



/***********************************************************************/
/* VERTEX ID MANAGEMENT */

	private:
		/* No. of vertex types in multigraph */
		static constexpr size_t n_vertex_types = std::tuple_size_v< vertex_types >;

		/* Array of free vertex IDs for automatic indexing */
		std::array<size_t, n_vertex_types> vertex_id_array;

		/* Register a vertex ID and increment ID counter */
		template<typename _VData>
		inline void register_vertex_id(const size_t id) {
			using check = assert_vertex_is_valid<_VData>;
			constexpr size_t I = vrdc::index_in_tuple< _VData, vertex_types>::value;
			assert(
				!has_vertex<_VData>(id) &&
				"Attempt to duplicate a vertex ID"
			);
			if (vertex_id_array[I] <= id) {
				vertex_id_array[I] = id + 1;
			}
		}

	protected:
		/* Get a new vertex ID */
		template<typename _VData>
		inline size_t next_vertex_id() {
			using check = assert_vertex_is_valid<_VData>;
			constexpr size_t I = vrdc::index_in_tuple<_VData, vertex_types>::value;
			return vertex_id_array[I];
		}

/***********************************************************************/
/* VERTEX LIST */

	private:
		/* Direct access to first and last vertices */
		template<typename _VData>
		struct as_vertex_list {
			using type = mlist< V<_VData> >;
		};

		using vertex_list_tuple_t = typename vrdc::transform< as_vertex_list, vertex_types >::type;

		vertex_list_tuple_t vertex_list_tuple;

	public:
		/* Access list of vertices */
		template<typename _VData>
		inline typename as_vertex_list<_VData>::type& vertices() {
			using check = assert_vertex_is_valid<_VData>;
			return std::get< typename as_vertex_list<_VData>::type >(vertex_list_tuple);
		}

/***********************************************************************/
/* VERTEX MAP */

	private:
		/* Tuple of vertex maps, one map per vertex type.
		Map key is the unique vertex index, map value is the vertex shared_ptr.
		Removal from map allows for vertex destruction once all other shared_ptr's
		are out of scope */
		template<typename _VData>
		struct as_vertex_map {
			using type = std::unordered_map< size_t, std::shared_ptr< V<_VData> > >;
		};

		using vertex_map_tuple_t = typename vrdc::transform< as_vertex_map, vertex_types >::type;

		vertex_map_tuple_t vertex_map_tuple;

		/* Access map of vertices (private, to avoid unsafe iteration) */
		template<typename _VData>
		inline typename as_vertex_map<_VData>::type& vertexMap() {
			using check = assert_vertex_is_valid<_VData>;
			return std::get< typename as_vertex_map<_VData>::type >(vertex_map_tuple);
		}

/***********************************************************************/
/* VERTEX LIFECYCLE */

	public:
		/* Insert vertex in graph, return raw pointer, known ID */
		template<typename _VData>
		inline V<_VData>* insert_vertex_ptr(const size_t id, _VData data = {}) {

			/* Keep track of provided vertex ID */
			register_vertex_id<_VData>(id);

			/* std::make_shared is not compatible with the
			protected constructor of V: get a raw ptr first */
			V<_VData>* vertex_ptr = new V<_VData>(*this, id, data);

			/* Register ID in vertex map */
			auto ret = vertexMap<_VData>().insert({ id, std::shared_ptr< V<_VData> >(vertex_ptr) });
			assert(ret.second == true);

			/* Append raw ptr to intrusive list */
			vertices<_VData>().push_back(vertex_ptr);

			return vertex_ptr;
		}

		/* Insert vertex in graph, return raw pointer, automatic ID */
		template<typename _VData>
		inline V<_VData>* insert_vertex_ptr(_VData data = {}) {
			const size_t id = next_vertex_id<_VData>();
			return insert_vertex_ptr<_VData>(id, data);
		}
	
		/* Insert vertex in graph, return shared_ptr, known ID */
		template<typename _VData>
		inline std::shared_ptr< V<_VData> > insert_vertex_shared_ptr(const size_t id, _VData data = {}) {

			V<_VData>* vertex_ptr = insert_vertex_ptr(id, data);

			return vertex_ptr->shared_from_this();
		}

		/* Insert vertex in graph, return shared_ptr, automatic ID*/
		template<typename _VData>
		inline std::shared_ptr < V<_VData> > insert_vertex_shared_ptr(_VData data = {}) {
			const size_t id = next_vertex_id<_VData>();
			return insert_vertex_shared_ptr<_VData>(id, data);
		}

		/* Insert vertex in graph, return reference, known ID */
		template<typename _VData>
		inline V<_VData>& insert_vertex(const size_t id, _VData data = {}) {
			return *insert_vertex_ptr(id, data);
		}

		/* Insert new vertex with automatic ID, return void */
		template<typename _VData>
		inline V<_VData>& insert_vertex(_VData data = {}) {
			const size_t id = next_vertex_id<_VData>();
			return insert_vertex<_VData>(id, data);
		}

		/* Delete a vertex and all its attached edges */
		template< typename _VData >
		inline void erase_vertex(V<_VData>& vertex) {
			/* Disconnect vertex: remove all edges and release shared_ptrs to the vertex */
			disconnect(vertex);
			
			if (vertex.active()) {
				/* Remove vertex from map: release last shared_ptr */
				size_t ret = vertexMap<_VData>().erase(vertex.id);
				assert(ret == 1 && "Attempt to delete an unregistered vertex.");
				/* Mark as inactive */
				vertex.active_ = false;
			}
		}

		/* Delete all edges to and from a vertex */
		template<typename _VData>
		void disconnect(V<_VData>& vertex) {
			/* Call function on unpacked tuple elements*/
			vrdc::foreach(vertex.edges_out_tuple, [&](auto& m) { clear_edge_list(m); });
			vrdc::foreach(vertex.edges_in_tuple, [&](auto& m) { clear_edge_list(m); });
		}

		/* Delete all vertices in a list */
		template<typename _List>
		void clear_vertex_list(_List& list) {
			static_assert(std::is_base_of_v<mlist_tag, _List>);
			for (auto& it : list) erase_vertex(it);
		}

		/* Delete all vertices in all lists */
		void reset() {
			vrdc::foreach(vertex_list_tuple, [&](auto& m) { clear_vertex_list(m); });
		}


/***********************************************************************/
/* VERTEX LOOKUP : FIND / HAS */

	public:
		/* Find vertex by ID, return shared_ptr, nullptr if not found */
		template<typename _VData>
		inline std::shared_ptr< V<_VData> > vertex_shared_ptr(const size_t id) {
			auto it = vertexMap<_VData>().find(id);
			if (it == vertexMap<_VData>().end()) {
				return nullptr;
			}
			return it->second;
		}

		/* Find vertex by ID, return raw ptr, nullptr if not found */
		template<typename _VData>
		inline V<_VData>* vertex_ptr(const size_t id) {
			auto it = vertexMap<_VData>().find(id);
			if (it == vertexMap<_VData>().end()) {
				return nullptr;
			}
			return it->second.get();
		}

		/* Find vertex by ID, return reference, throw out_of_range if not found */
		template<typename _VData>
		inline V<_VData>& vertex(const size_t id) {
			V<_VData>* ptr = vertex_ptr< _VData >(id);
			if (ptr == nullptr) {
				throw std::out_of_range("Vertex not found.");
			}
			return *ptr;
		}

		/* Check whether vertex ID exists */
		template<typename _VData>
		inline bool has_vertex(const size_t id) {
			return vertexMap<_VData>().find(id) != vertexMap<_VData>().end();
		}


/***********************************************************************/
/* EDGE MAP */

	private:
		/* No. of edge types in multigraph */
		static constexpr size_t n_edge_types = std::tuple_size_v< edge_types >;

		/* Tuple of edge maps, one map per edge type. Map key is the pair
		of indices of the source and target vertex */

		template<typename _Edge>
		using as_edge_map = std::unordered_map<
			edge_map_key const, std::shared_ptr< as_hgraph_edge< _Edge > > const,
			vrdc::hash_tuple< edge_map_key >
		>;

		template<typename... _Edge>
		using as_edge_map_tuple = std::tuple< as_edge_map< _Edge >... >;

		using edge_map_tuple_t = typename vrdc::repack_as_t<as_edge_map_tuple, edge_types>;
		edge_map_tuple_t edge_map_tuple;

		/* Acess map of edges (private, to avoid unsafe iteration) */
		template<typename _Edge>
		inline auto& edge_map() {
			using check = assert_edgetype_is_valid<_Edge>;
			return std::get< as_edge_map<_Edge> >(edge_map_tuple);
		}

/***********************************************************************/
/* EDGE LIFECYCLE */

		/* Remark ordering of template parameters: they follow the usual Edge ordering
		(source, target, data) when function arguments are IDs because the edge type must
		be specified as template, eg, myGraph::insert_edge< A, B, X >(0, 1). However,
		the data type parameter comes first when function arguments are references because
		source and target types are deduced, which allows for the more concise phrase
		myGraph::insert_edge< X >( myA, myB ). Data-less edges can omit the template specifier,
		eg, myGraph::insert_edge( myA, myB ).
		*/

	public:
		/* Insert edge in graph, return raw pointer */
		template< typename _EData = std::monostate, typename _Source, typename _Target >
		inline E<_Source, _Target, _EData>*
		insert_edge_ptr(V<_Source>& source, V<_Target>& target, _EData data = {}) {
			/* Shorthand for edge type */
			using edge_t = E<_Source, _Target, _EData>;
			/* Prepare map key */
			const auto edgeKey = std::make_tuple(source.id, target.id);
			/* Check uniqueness */
			assert(edge_map< edge_t >().find(edgeKey) == edge_map< edge_t >().end() &&
				"Attempt to duplicate an edge"
			);
			/* std::make_shared is not compatible with the
			protected constructor of V: get a raw ptr first */
			edge_t* edge_ptr = new edge_t(source, target, data);
			/* Register in map */
			auto ret = edge_map< edge_t >().insert({ edgeKey, std::shared_ptr< edge_t >(edge_ptr) });
			assert(ret.second == true);
			/* Append to intrusive lists in source and target */
			source.edges_to<_Target, _EData>().push_back(edge_ptr);
			target.edges_from<_Source, _EData>().push_back(edge_ptr);
			return edge_ptr;
		}

		/* Insert edge in graph, return raw pointer */
		template<typename _Source, typename _Target, typename _EData = std::monostate>
		inline E<_Source, _Target, _EData>*
		insert_edge_ptr(const size_t source_id, const size_t target_id, _EData data = {}) {
			/* Shorthand for edge type */
			using edge_t = E<_Source, _Target, _EData>;
			V<_Source>& source = vertex<_Source>(source_id);
			V<_Target>& target = vertex<_Target>(target_id);
			return insert_edge_ptr(source, target, data);
		}

		/* Insert edge in graph, return shared_ptr */
		template< typename _EData = std::monostate, typename _Source, typename _Target >
		inline std::shared_ptr< E<_Source, _Target, _EData> >
		insert_edge_shared_ptr(V<_Source>& source, V<_Target>& target,	_EData data = {}) {
			using edge_t = E<_Source, _Target, _EData>;
			edge_t* edge_ptr = insert_edge_ptr(source, target, data);
			return edge_ptr->shared_from_this();
		}

		/* Insert edge in graph, return shared_ptr */
		template<typename _Source, typename _Target, typename _EData = std::monostate>
		inline std::shared_ptr< E<_Source, _Target, _EData> >
		insert_edge_shared_ptr(const size_t source_id, const size_t target_id, _EData data = {}) {
			V<_Source>& source = vertex<_Source>(source_id);
			V<_Target>& target = vertex<_Target>(target_id);
			return insert_edge_shared_ptr(source, target, data);
		} 		

		/* Insert edge in graph, return reference */
		template<typename _EData = std::monostate, typename _Source, typename _Target>
		inline E<_Source, _Target, _EData>&
		insert_edge(V<_Source>& source, V<_Target>& target, _EData data = {}) {
			return *insert_edge_ptr(source, target, data);
		}

		/* Insert edge in graph, return reference */
		template<typename _Source, typename _Target, typename _EData = std::monostate>
		inline E<_Source, _Target, _EData>&
		insert_edge(const size_t source_id, const size_t target_id, _EData data = {}) {
			V<_Source>& source = vertex<_Source>(source_id);
			V<_Target>& target = vertex<_Target>(target_id);
			return insert_edge(source, target, data);
		}

		/* Delete an edge */
		template<typename _Source, typename _Target, typename _EData>
		inline void erase_edge(E<_Source, _Target, _EData>& edge) {
			using edge_t = E<_Source, _Target, _EData>;

			if (edge.active_) {
				/* Remove edge from graph's lookup map: releases shared_ptr */
				const size_t ret = edge_map< edge_t >().erase({ edge.source->id, edge.target->id });
				assert(ret == 1 && "Attempt to delete an unregistered edge.");
				/* Mark as inactive */
				edge.active_ = false;
			}
		}

		/* Delete all edges in an edge list */
		template<typename _List>
		void clear_edge_list(_List& list) {
			static_assert(std::is_base_of_v<mlist_tag, _List>);
			for (auto& it : list) erase_edge(it);
		}

/***********************************************************************/
/* EDGE LOOKUP : FIND / HAS */

		/* Remark ordering of template parameters: they follow the usual Edge ordering
		(source, target, data) when function arguments are IDs because the edge type must
		be specified as template, eg, myGraph::edge< A, B, X >(0, 1). However,
		the data type parameter comes first when function arguments are references because
		source and target types are deduced, which allows for the more concise phrase
		myGraph::edge< X >( myA, myB ). Data-less edges can omit the template specifier,
		eg, myGraph::edge( myA, myB ).
		*/

	public:
		/* Find edge by vertex IDs, return shared_ptr, nullptr if not found */
		template<typename _Source, typename _Target, typename _EData = std::monostate>
		inline std::shared_ptr< E<_Source, _Target, _EData> >
		edge_shared_ptr(const size_t source_id, const size_t target_id) {
			using edge_t = E<_Source, _Target, _EData>;
			auto it = edge_map< edge_t >().find({ source_id, target_id });
			if (it == edge_map< edge_t >().end()) {
				return nullptr;
			}
			return it->second;
		}

		/* Find edge by vertex references, return shared_ptr, nullptr if not found */
		template<typename _EData = std::monostate, typename _Source, typename _Target>
		inline std::shared_ptr< E<_Source, _Target, _EData> >
		edge_shared_ptr(const V<_Source>& source, const V<_Target>& target) {
			return edge_shared_ptr<_Source, _Target, _EData>(source.id, target.id);
		}

		/* Find edge by vertex IDs, return raw pointer, nullptr if not found */
		template<typename _Source, typename _Target, typename _EData = std::monostate>
		inline E<_Source, _Target, _EData>*
		edge_ptr(const size_t source_id, const size_t target_id) {
			using edge_t = E<_Source, _Target, _EData>;
			auto it = edge_map< edge_t >().find({ source_id, target_id });
			if (it == edge_map< edge_t >().end()) {
				return nullptr;
			}
			return it->second.get();
		}

		/* Find edge by vertex references, return raw pointer, nullptr if not found */
		template<typename _EData = std::monostate, typename _Source, typename _Target>
		inline E<_Source, _Target, _EData>*
		edge_ptr(const V<_Source>& source, const V<_Target>& target) {
			return edge_ptr<_Source, _Target, _EData>(source.id, target.id);
		}

		/* Find edge by vertex IDs, return reference, throw out_of_range if not found */
		template<typename _Source, typename _Target, typename _EData = std::monostate>
		inline E<_Source, _Target, _EData>&
		edge(const size_t source_id, const size_t target_id) {
			using edge_t = E<_Source, _Target, _EData>;
			auto it = edge_map< edge_t >().find({ source_id, target_id });
			if (it == edge_map< edge_t >().end()) {
				throw std::out_of_range("Edge not found.");
			}
			return *it->second;
		}

		/* Find edge by vertex references, return reference, throw out_of_range if not found */
		template<typename _EData = std::monostate, typename _Source, typename _Target>
		inline E<_Source, _Target, _EData>&
		edge(const V<_Source>& source, const V<_Target>& target) {
			return edge<_Source, _Target, _EData>(source.id, target.id);
		}

		/* Check whether edge exists, by vertex IDs */
		template<typename _Source, typename _Target, typename _EData = std::monostate>
		inline bool has_edge(const size_t source_id, const size_t target_id) {
			using edge_t = E<_Source, _Target, _EData>;
			return edge_map< edge_t >().find({ source_id, target_id }) != edge_map< edge_t >().end();
		}

		/* Check whether edge exists, by vertex reference */
		template<typename _EData = std::monostate, typename _Source, typename _Target>
		inline bool has_edge(const V<_Source>& source, const V<_Target>& target) {
			return has_edge<_Source, _Target, _EData>(source.id, target.id);
		}

/***********************************************************************/
/* TRAVERSAL */

	private:
		/***********************************************/
		/* Predicates for Variadic::foreach_if */

		/* Predicate: edge has data flag */
		template<typename _Edge, typename _EData>
		struct edge_has_type : std::is_base_of<_EData, typename _Edge::data_t> {};

		template<typename _IList, typename _EData>
		struct edgelist_has_type : edge_has_type<typename _IList::type, _EData> {};

	public:
		/* Traverse a tuple of edge maps and apply function for edges */
		template<typename _EData, typename _Lambda, typename... Ttup, typename... A>
		inline void foreach_edge_type(std::tuple<Ttup...>& tup, _Lambda lambda, A... args) {
			return vrdc::foreach_if<edgelist_has_type, _EData>(tup, lambda, args...);
		};

/***********************************************************************/
/* DIAGNOSTICS */

	public:
		
		
		/* No. of vertices in graph */
		size_t order() {
			size_t n_vertices = 0;
			vrdc::foreach(vertex_map_tuple, [&n_vertices](auto& map) {
				n_vertices += map.size();
				});
			return n_vertices;
		}

		/* Graph size as std::pair<no. of vertices, no. of edges> */
		std::pair<size_t, size_t> size() {
			size_t n_vertices = 0;
			size_t n_edges = 0;
			vrdc::foreach(vertex_map_tuple, [&n_vertices, &n_edges](auto& nodemap) {
				n_vertices += nodemap.size();
				for (auto it : nodemap) {
					vrdc::foreach(it.second->edges_out_tuple, [&n_edges](auto& edgemap) {
						n_edges += edgemap.size();
						});
				}
				});
			return { n_vertices, n_edges };
		}

/***********************************************************************/
/* RANDOM SHUFFLING */

	public:

		void shuffle() {

			/* Shuffle vertex lists */
			vrdc::foreach(vertex_list_tuple, [&](auto& m) { shuffle_list(m); });

			/* Shuffle edge lists */
			vrdc::foreach(vertex_list_tuple, [&](auto& m) { shuffle_edges(m); });
		}

		template<typename _List>
		void shuffle_list(_List& list) {
			static_assert(std::is_base_of_v<mlist_tag, _List>);

			/* Shuffle vector of pointers */
			auto rv = list.as_vector();			
			std::random_device rd;
			std::mt19937 g(rd());
			std::shuffle(rv.begin(), rv.end(), g);

			list.reorder(rv);
		}

		template<typename _List>
		void shuffle_edges(_List& list) {
			static_assert(std::is_base_of_v<mlist_tag, _List>);

			/* Scan vertices in list and shuffle all edges */
			for (auto& it : list) {
				vrdc::foreach(it.edges_out_tuple, [&](auto& m) { shuffle_list(m); });
				vrdc::foreach(it.edges_in_tuple, [&](auto& m) { shuffle_list(m); });
			}
		}

/***********************************************************************/
/* CONSTRUCTOR / DESTRUCTOR */

	public:
		graph() : vertex_id_array{} {}
	};
}