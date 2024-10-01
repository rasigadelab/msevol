#pragma once
#include "../hgraph/variadic.h"

namespace hgraph{
	template<typename _Graph>
	struct meta {
		/* Vertex */
		template<typename _VData>
		using V = typename _Graph::template V< _VData >;

		/* Edge */
		template<typename _Source, typename _Target, typename _EData>
		using E = typename _Graph::template  E< _Source, _Target, _EData >;

		/* Declared vertex types */
		using vertex_types = typename _Graph::vertex_types;

		/* Declared edge types */
		using edge_types = typename _Graph::edge_types;

		/* Edge type manipulation*/
		template<typename _Source, typename _Target>
		using as_target = _Target;

		template<typename _Source, typename _Target>
		using as_source = _Source;

		/* Transform a tuple of edge types */
		template<template<typename, typename> typename, typename...> struct transform_edges;
		template<template<typename, typename> typename T, typename... Ts>
		struct transform_edges < T, std::tuple < Ts... >> {
			using type = std::tuple< T< typename Ts::source_t, typename Ts::target_t>... >;
		};

		/* Multiplicity edge types */
		template<typename _Edge>
		struct is_multiplicity_edge : std::is_base_of<Mlt, typename _Edge::data_t> {};

		/* Tuple type of multiplicity edges */
		using multiplicity_edge_types = typename vrdc::repack_if<is_multiplicity_edge, edge_types>::type;

		/* Containers of a vertex type */
		template<typename _VData>
		struct container_types {
			using multiplicity_edges_out_t = typename vrdc::repack_if< is_multiplicity_edge, typename V<_VData>::edges_out_tuple_t >::type;
			using type = typename transform_edges< as_target, multiplicity_edges_out_t >::type;
		};

		/* Containees of a vertex type */
		template<typename _VData>
		struct containee_types {
			using multiplicity_edges_in_t = typename vrdc::repack_if< is_multiplicity_edge, typename V<_VData>::edges_in_tuple_t >::type;
			using type = typename transform_edges< as_source, multiplicity_edges_in_t >::type;
		};
	};
}
