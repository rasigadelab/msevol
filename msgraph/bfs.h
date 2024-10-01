#pragma once

/*********************************************************************************/
/*	BREADTH-FIRST TRAVERSAL OF MSGRAPH

See bfs_callback for an exemplary callback structure .

bfs_queue is used internally by BFS class.

BFS class implements traversal along multiplicity edges.

*/
/*********************************************************************************/
#include <queue>
#include <unordered_set>

#include "meta.h"
#include "../hgraph/variadic.h"

namespace hgraph {

	/* BFS callback, exposes visit_vertex and visit_edge members */
	struct bfs_callback {
		template<typename T>
		void visit_first(T& v) {};

		template<typename T>
		void visit_vertex(T& v) {};

		template<typename T>
		void visit_edge_out(T& e) {};

		template<typename T>
		void visit_edge_in(T& e) {};
	};

	/**********************************************************************************/
	/* Queue structure for BFS traversal */
	template<typename _Graph>
	class bfs_queue {
	public:
		/* Template inheritance shorthands */
		template<typename _VData>
		using V = typename _Graph::template V< _VData >;
		using vertex_types = typename _Graph::vertex_types;

		/* Append to queue and mark as visited */
		template<typename _VData>
		inline void enqueue(V< _VData >& v) {
			queue< _VData >().push(&v);
			discovered_vertices< _VData >().insert(&v);
		}

		/* Pop first vertex */
		template<typename _VData>
		inline V< _VData >* dequeue() {
			auto& q = queue< _VData >();
			V< _VData >* vertex_ptr = q.front();
			q.pop();
			return vertex_ptr;
		}

		/* Check if vertex has been visited */
		template<typename _VData>
		inline bool is_discovered(V< _VData >& v) {
			auto& s = discovered_vertices< _VData >();
			return s.find(&v) != s.end();
		}

		/* Queue size */
		template<typename _VData>
		inline size_t size() {
			return queue< _VData >().size();
		}

	protected:
		/* Queue of vertices to visit */
		template<typename _VData>
		using queue_t = std::queue< V< _VData >* >;

		using queue_tuple_t = typename vrdc::transform_as< queue_t, vertex_types >::type;
		queue_tuple_t queue_tuple;

		template<typename _VData>
		inline queue_t< _VData >& queue() {
			return std::get< queue_t< _VData > >(queue_tuple);
		}

		/* Discovered vertices */
		template<typename _VData>
		using set_t = std::unordered_set< V< _VData >* >;

		using set_tuple_t = typename vrdc::transform_as< set_t, vertex_types >::type;
		set_tuple_t set_tuple;

		template<typename _VData>
		inline set_t< _VData >& discovered_vertices() {
			return std::get< set_t< _VData > >(set_tuple);
		}

		template<typename _VData>
		inline void mark_as_discovered(V< _VData >& v) {
			discovered_vertices< _VData >().insert(&v);
		}
	};

	/* BFS direction options: scan edges forward or backward ? */
	enum class bfs_direction { out,	in };

	/**********************************************************************************/
	/* Implements breadth-first traversal of multiplicity graph. */
	template<typename _Graph, bfs_direction _Direction = bfs_direction::out>
	class BFS {
	public:
		/* Template inheritance shorthands */
		template<typename _VData>
		using V = typename _Graph::template V< _VData >;
		template<typename _Source, typename _Target, typename _EData>
		using E = typename _Graph::template E< _Source, _Target, _EData >;
		using vertex_types = typename _Graph::vertex_types;
		using edge_types = typename _Graph::edge_types;

		/* Interface */
		BFS(_Graph& graph) : graph{ graph }, Q{ bfs_queue< _Graph >() } {};

		template<typename _VData, typename _Callback = std::monostate>
		void run(V< _VData >& v, _Callback& cb = {}) {
			/* Enqueue starting vertex and start scan */
			Q = bfs_queue< _Graph >();
			Q.enqueue(v);
			/* Launch callback */
			if constexpr (!std::is_same_v< _Callback, std::monostate >) {
				cb.visit_first(v);
			}
			static constexpr size_t Istart = vrdc::index_in_tuple< _VData, vertex_types >::value;
			scan_vertices< Istart >(cb);
		}

	protected:
		_Graph& graph;
		bfs_queue< _Graph > Q;

		/* Scan outgoing edges */
		template<size_t I = 0, typename _VData, typename _Callback>
		void scan_edges_out(V< _VData >& v, _Callback& cb) {

			using container_types = typename meta< _Graph >::template container_types< _VData >::type;

			if constexpr (I < std::tuple_size_v< container_types >) {
				using container_t = std::tuple_element_t< I, container_types >;

				/* Directional switch */
				for (auto& e : v.edges_to< container_t, Mlt >()) {

					/* Launch callback */
					if constexpr (!std::is_same_v< _Callback, std::monostate >) {
						cb.visit_edge_out(e);
					}

					auto& container = *e.target;
					if (!Q.is_discovered(container)) {

						/* Enqueue container and mark as visited */
						Q.enqueue(container);
					}

				}
				scan_edges_out< I + 1 >(v, cb);
			}
		}

		/* Scan incoming edges */
		template<size_t I = 0, typename _VData, typename _Callback>
		void scan_edges_in(V< _VData >& v, _Callback& cb) {

			using containee_types = typename meta< _Graph >::template containee_types< _VData >::type;

			if constexpr (I < std::tuple_size_v< containee_types >) {
				using containee_t = std::tuple_element_t< I, containee_types >;

				/* Directional switch */
				for (auto& e : v.edges_from< containee_t, Mlt >()) {

					/* Launch callback */
					if constexpr (!std::is_same_v< _Callback, std::monostate >) {
						cb.visit_edge_in(e);
					}

					auto& containee = *e.source;
					if (!Q.is_discovered(containee)) {

						/* Enqueue containee and mark as visited */
						Q.enqueue(containee);
					}

				}
				scan_edges_in< I + 1 >(v, cb);
			}
		}

		/* Iterate over toposorted vertex types. */
		template<size_t I, typename _Callback>
		void scan_vertices(_Callback& cb) {

			if constexpr (I < std::tuple_size_v< vertex_types >) {
				using data_t = std::tuple_element_t< I, vertex_types >;

				while (Q.size< data_t >() > 0) {

					/* Dequeue vertex */
					V< data_t >* vertex_ptr = Q.dequeue< data_t >();

					/* Launch callback */
					if constexpr (!std::is_same_v< _Callback, std::monostate >) {
						cb.visit_vertex(*vertex_ptr);
					}
					
					if constexpr (_Direction == bfs_direction::out) {
						scan_edges_out(*vertex_ptr, cb);
					}
					else {
						scan_edges_in(*vertex_ptr, cb);
					}
				}

				if constexpr (_Direction == bfs_direction::out) {
					scan_vertices< I + 1 >(cb);
				}
				else {
					scan_vertices< I - 1 >(cb);
				}
			}
		}
	};	
}