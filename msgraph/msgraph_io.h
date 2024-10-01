#pragma once

/*******************************************************************/
/*
Serialization routines based on CSV files.
	- specify reflection/formatting info using _Reflect template parameter
*/
/*******************************************************************/

#include "../msgraph/io/Writer.h"
#include "../msgraph/io/Reader.h"
#include "../msgraph/msgraph.h"

template< typename _Reflect >
class hgraph_io {
public:

	/* Simplify template argument calls */
	template< typename T >
	using reflect = typename _Reflect::template get<T>;

	/* Graph type */
	using graph_type = typename _Reflect::graph_type;
	using state_type = typename _Reflect::ModelState;

	using state_get = typename _Reflect::template get< state_type >;

	/* Vertex shorthand */
	template< typename _VData >
	using V = typename graph_type::template V< _VData >;

	/* Edge shorthand */
	template<typename _Source, typename _Target, typename _EData>
	using E = typename graph_type::template E<_Source, _Target, _EData>;

	/*******************************************************************/
	/* VERTEX SERIALIZATION: WRITERS */

		/* Generic vertex serialization info */
	template< typename _VData >
	struct getVertex {
		static std::string name() {
			std::string vertex_name(reflect<_VData>::name);
			return std::string("Vertex_" + vertex_name);
		}
		static auto header_tuple() {
			return std::tuple_cat(
				std::make_tuple("id"),
				reflect< _VData >::header_tuple()
			);
		}
		static auto as_tuple(V< _VData >& x) {
			return std::tuple_cat(
				std::make_tuple(x.id),
				reflect< _VData >::as_tuple(x.data)
			);
		}
	};

	/* Write a vertex */
	template< typename _VData > class VertexWriter : public Writer {
	public:
		using type = _VData;
		using autoname_t = getVertex<_VData>;
		VertexWriter(const char* dir, const bool append = false) :
			Writer(dir, autoname_t::name().c_str(), append) {
			if (!append) {
				firstline_tuple(std::tuple_cat(autoname_t::header_tuple(), state_get::header_tuple()));
			}
		}
		void put(V< _VData >& x, const state_type state = {}) {
			push_tuple(std::tuple_cat( autoname_t::as_tuple(x), state_get::as_tuple(state) ));
		}
	};

	/* Write all vertices of given type */
	template<typename _VData>
	static void writeVertices(graph_type& u, const char* dataDir, const state_type state = {}, const bool append = false) {
		VertexWriter< _VData > writer(dataDir, append);
		for (auto& vertex : u.vertices< _VData >()) {
			writer.put(vertex, state);
		}
	}

	/* Write all vertices of all types */
	template<size_t I = 0>
	static void writeAllVertices(graph_type& u, const char* dataDir, const state_type state = {}, const bool append = false) {
		using type = std::tuple_element_t < I, typename graph_type::vertex_types > ;
		static constexpr bool is_serializable = !std::is_same_v< typename reflect< type >::input_tuple_t, std::false_type >;
		if constexpr (is_serializable) {
			writeVertices< type >(u, dataDir, state, append);
		}
		if constexpr (I + 1 < std::tuple_size_v< typename graph_type::vertex_types >) {
			writeAllVertices< I + 1 >(u, dataDir, state, append);
		}
	}

	/*******************************************************************/
	/* VERTEX SERIALIZATION: READERS */

	/* Read a vertex */
	template< typename _VData > class VertexReader : public Reader {
	public:
		using type = _VData;
		using autoname_t = getVertex<_VData>;
		using tuple_t = typename reflect<_VData>::input_tuple_t;

		VertexReader(const char* dir) :
			Reader(dir, autoname_t::name().c_str()) {}

		bool fetch(graph_type& u) {
			size_t id = 0;
			tuple_t t;
			if (read_tuple(id, t)) {
				_VData data = std::make_from_tuple< _VData >(std::move(t));
				auto vertex = u.insert_vertex_shared_ptr(id, data);
				return true;
			}
			else return false;
		}
	};

	/* Read all vertices of given type */
	template<typename _VData>
	static void readVertices(graph_type& u, const char* dataDir) {
		VertexReader< _VData > reader(dataDir);
		while (reader.fetch(u)) {}
	}

	/* Read all vertices of all types */
	template<size_t I = 0>
	static void readAllVertices(graph_type& u, const char* dataDir) {
		using type = std::tuple_element_t< I, typename graph_type::vertex_types >;
		static constexpr bool is_serializable = !std::is_same_v< typename reflect< type >::input_tuple_t, std::false_type >;
		if constexpr (is_serializable) {
			readVertices< type >(u, dataDir);
		}
		if constexpr (I + 1 < std::tuple_size_v< typename graph_type::vertex_types >) {
			readAllVertices< I + 1 >(u, dataDir);
		}
	}

	/*******************************************************************/
	/* EDGE SERIALIZATION: WRITERS */

	/* Generic edge serialization info */
	template<typename _Source, typename _Target, typename _EData>
	struct getEdge {
		static std::string name() {
			std::string source_name(reflect< _Source >::name);
			std::string target_name(reflect< _Target >::name);
			std::string data_name(reflect< _EData >::name);
			std::string edge_name;
			edge_name += "Edge_" + source_name + "_" + target_name + "_" + data_name;
			return edge_name;
		}
		static auto header_tuple() {
			return std::tuple_cat(
				std::make_tuple("source", "target"),
				reflect< _EData >::header_tuple()
			);
		}
		static auto as_tuple(E< _Source, _Target, _EData >& x) {
			return std::tuple_cat(
				std::make_tuple(x.source->id, x.target->id),
				reflect< _EData >::as_tuple(x.data)
			);
		}
	};

	/* Write an edge */
	template<typename _Source, typename _Target, typename _EData>
	class EdgeWriter : public Writer {
	public:
		using source_t = _Source;
		using target_t = _Target;
		using data_t = _EData;
		using autoname_t = getEdge< source_t, target_t, data_t >;

		EdgeWriter(const char* dir, const bool append = false) :
			Writer(dir, autoname_t::name().c_str(), append) {
			if (!append) {
				firstline_tuple(std::tuple_cat(autoname_t::header_tuple(), state_get::header_tuple()));
			}
		}

		void put(E< _Source, _Target, _EData >& x, const state_type state = {}) {
			push_tuple(std::tuple_cat(autoname_t::as_tuple(x), state_get::as_tuple(state)));
		}
	};

	/* Write all edges of a given type */
	template<typename _Source, typename _Target, typename _Data>
	static void writeEdges(graph_type& u, const char* dataDir, const state_type state = {}, const bool append = false) {
		EdgeWriter< _Source, _Target, _Data > writer(dataDir, append);

		/* Scan each vertex */
		for (auto& vertex : u.vertices< _Source >()) {
			for (auto& edge : vertex.edges_to<_Target, _Data>()) {
				writer.put(edge, state);
			}
		}
	}

	/* Write all edges of all types */
	template<size_t I = 0>
	static void writeAllEdges(graph_type& u, const char* dataDir, const state_type state = {}, const bool append = false) {
		using edge_t = std::tuple_element_t< I, typename graph_type::edge_types >;
		using source_t = typename edge_t::source_t;
		using target_t = typename edge_t::target_t;
		using data_t = typename edge_t::data_t;

		static constexpr bool is_serializable = !std::is_same_v< typename reflect< data_t >::input_tuple_t, std::false_type >;
		if constexpr (is_serializable) {
			writeEdges< source_t, target_t, data_t >(u, dataDir, state, append);
		}
		if constexpr (I + 1 < std::tuple_size_v< typename graph_type::edge_types >) {
			writeAllEdges< I + 1 >(u, dataDir, state, append);
		}
	}

	/*******************************************************************/
	/* EDGE SERIALIZATION: READERS */

	/* Generic edge reader */
	template<typename _Source, typename _Target, typename _EData>
	class EdgeReader : public Reader {
	public:
		using source_t = _Source;
		using target_t = _Target;
		using data_t = _EData;
		using tuple_t = typename reflect<_EData>::input_tuple_t;
		using autoname_t = getEdge< source_t, target_t, data_t >;

		EdgeReader(const char* dir) :
			Reader(dir, autoname_t::name().c_str()) {}

		bool fetch(graph_type& u) {
			size_t source_id = 0, target_id = 0;
			tuple_t t;
			if (read_tuple(source_id, target_id, t)) {
				auto source = u.vertex_shared_ptr< _Source >(source_id);
				auto target = u.vertex_shared_ptr< _Target >(target_id);
				_EData data = std::make_from_tuple< _EData >(std::move(t));
				auto edge = u.insert_edge_shared_ptr(*source, *target, data);
				return true;
			}
			else return false;
		}
	};

	/* Specialization: multiplicity edge */
	template<typename _Source, typename _Target>
	class EdgeReader< _Source, _Target, Mlt > : public Reader {
	public:
		using source_t = _Source;
		using target_t = _Target;
		using data_t = Mlt;
		using tuple_t = typename reflect<Mlt>::input_tuple_t;
		using autoname_t = getEdge< source_t, target_t, data_t >;

		EdgeReader(const char* dir) :
			Reader(dir, autoname_t::name().c_str()) {}

		bool fetch(graph_type& u) {
			size_t source_id = 0, target_id = 0;
			tuple_t t;
			if (read_tuple(source_id, target_id, t)) {
				auto source = u.vertex_shared_ptr< _Source >(source_id);
				auto target = u.vertex_shared_ptr< _Target >(target_id);
				size_t data = std::make_from_tuple< size_t >(std::move(t));
				auto edge = u.assign_shared_ptr(*source, *target, data);
				return true;
			}
			else return false;
		}
	};

	/* Read all edges of a given type */
	template<typename _Source, typename _Target, typename _Data>
	static void readEdges(graph_type& u, const char* dataDir) {
		EdgeReader< _Source, _Target, _Data > reader(dataDir);
		while (reader.fetch(u)) {}
	}

	/* Read all edges of a all types */
	template<size_t I = 0>
	static void readAllEdges(graph_type& u, const char* dataDir) {
		using edge_t = std::tuple_element_t< I, typename graph_type::edge_types >;
		using source_t = typename edge_t::source_t;
		using target_t = typename edge_t::target_t;
		using data_t = typename edge_t::data_t;

		static constexpr bool is_serializable = !std::is_same_v< typename reflect< data_t >::input_tuple_t, std::false_type >;
		if constexpr (is_serializable) {
			readEdges< source_t, target_t, data_t >(u, dataDir);
		}
		if constexpr (I + 1 < std::tuple_size_v< typename graph_type::edge_types >) {
			readAllEdges< I + 1 >(u, dataDir);
		}
	}

	/*******************************************************************/
	/* FULL-GRAPH WRITE / READ */

	static void readGraph(graph_type& u, const char* dataDir) {
		readAllVertices(u, dataDir);
		readAllEdges(u, dataDir);
	}

	/*******************************************************************/
	/* FULL-GRAPH WRITE with step no. */

	static void writeGraph(graph_type& u, const char* dataDir, const state_type state = {}, const bool append = false) {
		writeAllVertices(u, dataDir, state, append);
		writeAllEdges(u, dataDir, state, append);
	}

}; // class hgraph_io