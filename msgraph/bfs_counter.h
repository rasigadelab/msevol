#pragma once

/*********************************************************************************/
/*	MULTISET COUNTER

Provide this as a callback to BFS traversal routine.

Upward count: traverses multiplicity edges from containees to containers. Stores
the count of the starting vertex (containee) in each visited container. Access stored
counts using ::count() method.

*/
/*********************************************************************************/

#include <unordered_map>
#include <vector>
#include <cassert>

namespace hgraph {

	template<typename _Graph>
	class bfs_toposort {
	protected:
		_Graph& graph;
	public:
		explicit bfs_toposort(_Graph& graph) : graph{ graph } {}

		/* Template inheritance shorthands */
		template<typename _VData>
		using V = typename _Graph::template V< _VData >;
		template<typename _Source, typename _Target, typename _EData>
		using E = typename _Graph::template E< _Source, _Target, _EData >;
		using vertex_types = typename _Graph::vertex_types;
		using edge_types = typename _Graph::edge_types;

		template< typename _VData >
		using vec_t = std::vector< V<_VData>* >;

		using vec_tuple_t = typename vrdc::transform_as< vec_t, vertex_types >::type;
		vec_tuple_t vec_tuple;

		template< typename _VData >
		inline vec_t< _VData >& vec() {
			return std::get< vec_t< _VData > >(vec_tuple);
		}

		template<typename _VData>
		void visit_first(V<_VData>& v) {}

		template<typename _VData>
		void visit_vertex(V<_VData>& v) {
			vec< _VData >().push_back(&v);
		}

		template<typename _Source, typename _Target>
		void visit_edge_out(E<_Source, _Target, Mlt>& e) {}

		template<typename _Source, typename _Target>
		void visit_edge_in(E<_Source, _Target, Mlt>& e) {}

	};



	template<typename _Graph>
	class bfs_counter {
	protected:
		_Graph& graph; 
	public:
		bfs_counter(_Graph& graph) : graph{ graph } {}


			/* Tuple functor: copy */
		template<typename T>
		struct tup_copy {
			inline void operator()(T& target_map, bfs_counter< _Graph >& source) {
				target_map = std::get< T >(source.counter_map_tuple);
			}
		};

		bfs_counter< _Graph >& operator=(bfs_counter< _Graph >& copy) {
			graph = copy.graph;
			vrdc::tup_apply< tup_copy >(counter_map_tuple, copy);
			return *this;
		}

		bfs_counter(bfs_counter< _Graph >& copy) : graph{ copy.graph } {
			vrdc::tup_apply< tup_copy >(counter_map_tuple, copy);
		}

		/* Template inheritance shorthands */
		template<typename _VData>
		using V = typename _Graph::template V< _VData >;
		template<typename _Source, typename _Target, typename _EData>
		using E = typename _Graph::template E< _Source, _Target, _EData >;
		using vertex_types = typename _Graph::vertex_types;
		using edge_types = typename _Graph::edge_types;

		template< typename _VData >
		using counter_map = std::unordered_map< V<_VData>*, size_t >;

		using counter_map_tuple_t = typename vrdc::transform_as< counter_map, vertex_types >::type;
		counter_map_tuple_t counter_map_tuple;

		template< typename _VData >
		inline counter_map< _VData >& map() {
			return std::get< counter_map< _VData > >(counter_map_tuple);
		}

		template< typename _VData >
		inline const counter_map< _VData > cmap() {
			return std::get< counter_map< _VData > >(counter_map_tuple);
		}

		/* Callbacks */
		template<typename _VData>
		void visit_first(V<_VData>& v) {
			auto& cmap = map< _VData >();
			cmap.insert({ &v, 1 });
		}

		template<typename _VData>
		void visit_vertex(V<_VData>& v) {}

		template<typename _Source, typename _Target>
		void visit_edge_out(E<_Source, _Target, Mlt>& e) {
			/* At each edge visit, update count of target */
			auto& cmap_source = map< _Source >();
			auto& cmap_target = map< _Target >();

			auto it_source = cmap_source.find(e.source.get());
			assert(it_source != cmap_source.end());

			size_t count = it_source->second * e.data.mult;

			auto it_target = cmap_target.find(e.target.get());
			if (it_target == cmap_target.end()) {
				auto ret = cmap_target.insert({ e.target.get(), 0 });
				assert(ret.second);
				it_target = ret.first;
			}
			it_target->second += count;
		}

		template<typename _Source, typename _Target>
		void visit_edge_in(E<_Source, _Target, Mlt>& e) {
			/* At each edge visit, update count of source */
			auto& cmap_source = map< _Source >();
			auto& cmap_target = map< _Target >();

			auto it_target = cmap_target.find(e.target.get());
			assert(it_target != cmap_target.end());

			size_t count = it_target->second * e.data.mult;

			auto it_source = cmap_source.find(e.source.get());
			if (it_source == cmap_source.end()) {
				auto ret = cmap_source.insert({ e.source.get(), 0 });
				assert(ret.second);
				it_source = ret.first;
			}
			it_source->second += count;
		}

		/* Return count associated with the given vertex; 0 if not found */
		template<typename _VData>
		size_t count(V<_VData>& v) {
			auto& cmap = map< _VData >();
			auto it = cmap.find(&v);
			if (it == cmap.end()) {
				return 0;
			}
			else {
				return it->second;
			}
		}

		/* Starting vertex */
		//V<_VData>* start;


	};
}
