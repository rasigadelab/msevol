#pragma once

/******************************************************************/
/* MULTISET */

/* Generic structures to store multiplicity changes and barcodes */

#include <tuple>
#include <algorithm>

#include "meta.h"

namespace hgraph {

	template<typename _Graph>
	class multiset {
	public:
		/* Template inheritance shorthands */
		template<typename _VData>
		using V = typename _Graph::template V< _VData >;
		template<typename _Source, typename _Target, typename _EData>
		using E = typename _Graph::template E< _Source, _Target, _EData >;
		using vertex_types = typename _Graph::vertex_types;
		using edge_types = typename _Graph::edge_types;

		/* Find index of vertex type in the list of types assigned to the graph. Used to initialize
			hash values depending to vertex types. */
		template<typename _VData>
		static constexpr size_t vertex_type_index() {
			static constexpr size_t I = vrdc::index_in_tuple< _VData, vertex_types >::value;
			static_assert(I < std::tuple_size_v< vertex_types >, "Vertex type not declared in graph.");
			return I;
		}
			   
		/* Multiset element holding pointer to containee and multiplicity. Pointer can be nullptr
		to represent destruction of an element in dispatching algorithms. */
		template< typename _VData >
		struct mset_element {
			V< _VData >* vertex_ptr;
			size_t mult;
			mset_element(V< _VData >* const vertex_ptr, size_t mult) :
				vertex_ptr{ vertex_ptr }, mult{ mult } {}

			/* Lexicographic ordering: y >= x ignoring multiplicity */
			static inline bool notless_by_id(const mset_element<_VData>& x, const mset_element<_VData>& y) {
				return x.vertex_ptr < y.vertex_ptr;
			}

			inline bool operator==(const mset_element< _VData >& x) const {
				return (vertex_ptr == x.vertex_ptr) && (mult == x.mult);
			}

			inline bool operator!=(const mset_element< _VData >& x) const {
				return !(x == *this);
			}
		};

		/* Homogeneous multiset component representing a multiset of containees of the same type.
		Combine these to obtain an heterogeneous multiset. Containee uniqueness is not enforced upon
		insertion, responsibility for uniqueness is left to the caller. */
		template< typename _VData >
		struct mset_vector : std::vector< mset_element< _VData > > {
			using data_t = _VData;

			/* Extract vector of multiplicities for lower-level routines */
			std::vector< size_t > get_mult_vector() {
				std::vector< size_t > v(this->size());
				for (size_t i = 0; i < this->size(); ++i) {
					v[i] = (*this)[i].mult;
				}
				return v;
			}
			inline mset_vector< _VData >& operator-=(const std::vector< size_t >& v) {
				assert(this->size() == v.size());
				for (size_t i = 0; i < this->size(); ++i) {
					(*this)[i].mult -= v[i];
				}
				return *this;
			}
			/* Replace multiplicities with a those of a vector. */
			inline mset_vector< _VData >& set_mult_from(const std::vector< size_t >& v) {
				assert(this->size() == v.size());
				for (size_t i = 0; i < this->size(); ++i) {
					(*this)[i].mult = v[i];
				}
				return *this;
			}
			/* Sort by increasing pointer values. */
			inline void sort() {
				std::sort(this->begin(), this->end(), mset_element< _VData >::notless_by_id);
			}
			/* Hash function called by mset's hash() functor. Quick'n'dirty implementation
			with much room for improvement. */
			inline size_t hash() const {
				/* Hash initialized depending on containee type to avoid hash collision between
				empty multisets of different types */
				size_t h = vertex_type_index< _VData >() << 7;
				for (size_t i = 0; i < this->size(); ++i) {
					h ^= std::hash< size_t >()(size_t((*this)[i].vertex_ptr) << 1) ^
						std::hash< size_t >()((*this)[i].mult << 3);
				}
				return h;
			}
		};

		/* Heterogeneous multiset representation. Only defined for population containers. The multiset
		is represented as a tuple of multiset vectors with one vector per type of containee. Use vec<T>()
		to access such vectors. */
		template<typename _VData>
		struct mset {
			static_assert(std::is_base_of_v<population_container, _VData>);


			/* Repack a multiplicity edge as mset_vector */
			template<typename _Source, typename _Target>
			using as_mset_vector = mset_vector< _Source >;

			template<typename _Source, typename _Target>
			using as_source = _Source;

			/* A set of vectors of possible containees of Vdata */
			using multiplicity_edges_in_t = typename vrdc::repack_if< typename meta<_Graph>::is_multiplicity_edge, typename V<_VData>::edges_in_tuple_t >::type;
			using containee_types = typename meta< _Graph>::template transform_edges< as_source, multiplicity_edges_in_t >::type;
			using mset_vector_tuple_t = typename vrdc::transform_as< mset_vector, containee_types >::type;
			mset_vector_tuple_t mset_vector_tuple;

			/* Non-const access to homogeneous multiset vector for a containee type */
			template<typename _Source>
			inline mset_vector< _Source >& vec() {
				return std::get< mset_vector< _Source > >(mset_vector_tuple);
			}
			/* Const access to homogeneous multiset vector for a containee type */
			template<typename _Source>
			inline const mset_vector< _Source >& cvec() const {
				return std::get< mset_vector< _Source > >(mset_vector_tuple);
			}
			/* Scan containees of a multiset vertex to fill the multiset vectors. Resets all vectors before scanning. */
			template<typename _Source>
			inline void scan_specific(V<_VData>& v) {
				auto& edges = v.edges_from< _Source, Mlt >();
				auto& vect = vec< _Source >();
				if (vect.size() > 0) {
					vect.clear();
				}
				vect.reserve(edges.size());
				for (auto& edge : edges) {
					vect.emplace_back(edge.source.get(), edge.data.mult);
				}
				vect.sort();
			}
			/* Scan containees of a multiset vertex to fill the multiset vectors. Resets all vectors before scanning. */
			template<size_t I = 0>
			inline void scan(V<_VData>& v) {
				if constexpr (I < std::tuple_size_v<mset_vector_tuple_t>) {
					using data_t = typename std::tuple_element_t< I, mset_vector_tuple_t >::data_t;
					scan_specific< data_t >(v);
					scan< I + 1 >(v);
				}
			}

			/* Hash functor for use with unordered_map. */
			struct hash {
				template<typename _Source>
				static inline void hash_specific(const mset< _VData >& m, size_t& h) {
					h ^= m.cvec< _Source >().hash();
				}
				template<size_t I = 0>
				static inline void hash_generic(const mset< _VData >& m, size_t& h) {
					if constexpr (I < std::tuple_size_v<mset_vector_tuple_t>) {
						using data_t = typename std::tuple_element_t< I, mset_vector_tuple_t >::data_t;
						hash_specific< data_t >(m, h);
						hash_generic< I + 1 >(m, h);
					}
				}
				size_t operator()(const mset< _VData >& m) const {
					size_t h = vertex_type_index< _VData >() << 7;
					hash_generic(m, h);
					return h;
				}
			};
			/* Equality operators */
			template<size_t I = 0>
			inline void equal_generic(const mset< _VData >& x, bool& eq) const {
				if constexpr (I < std::tuple_size_v<mset_vector_tuple_t>) {
					if (eq) {
						using data_t = typename std::tuple_element_t< I, mset_vector_tuple_t >::data_t;
						const auto& lvec = this->cvec< data_t >();
						const auto& rvec = x.cvec< data_t >();
						eq = eq && (lvec == rvec);
						equal_generic< I + 1 >(x, eq);
					}
				}
			}
			inline bool operator==(const mset< _VData >& x) const {
				bool eq = true;
				equal_generic(x, eq);
				return eq;
			}
			inline bool operator!=(const mset< _VData >& x) const {
				return !(x == *this);
			}

			/* Constructors */
			mset(V<_VData>& v) {
				scan(v);
			}
			mset() {}

			/* Change multiplicity of an element. Inserts vector element if required, and erases
			an existing element if mult = 0. */
			template<typename _Source>
			void update(V< _Source >& source, const size_t mult) {
				static_assert(vrdc::has_type< _Source, containee_types >::value);

				auto& vect = vec< _Source >();
				mset_element< _Source > me(&source, mult);

				auto it = std::lower_bound(vect.begin(), vect.end(), me, mset_element< _Source >::notless_by_id);

				if (it != vect.end() && !(me.vertex_ptr < it->vertex_ptr)) {

					if (mult == 0) {
						/* Delete entry */
						vect.erase(it);
					}
					else {
						/* Update multiplicity */
						it->mult = mult;
					}
				}
				else {
					/* Insert sorted */
					if (mult > 0) {
						vect.insert(it, std::move(me));
					}
				}
			}

			/* Update several elements by merging with an outcome vector, ignoring nullptr entries. */
			template<typename _Source>
			void update(mset_vector< _Source >& v) {
				for (auto& it : v) {
					if (it.vertex_ptr != nullptr) {
						update(*it.vertex_ptr, it.mult);
					}
				}
			}

			/* Shift multiplicity of an element */
			template<typename _Source>
			void shift_mult(V< _Source >& source, const int64_t delta) {
				static_assert(vrdc::has_type< _Source, containee_types >::value);

				if (delta == 0) return;

				auto& vect = vec< _Source >();
				mset_element< _Source > me(&source, size_t(delta));

				auto it = std::lower_bound(vect.begin(), vect.end(), me, mset_element< _Source >::notless_by_id);

				if (it != vect.end() && !(me.vertex_ptr < it->vertex_ptr)) {
					assert(size_t(delta) <= it->mult);
					const size_t mult = it->mult + delta;

					if (mult == 0) {
						/* Delete entry */
						vect.erase(it);
					}
					else {
						/* Update multiplicity */
						it->mult = mult;
					}
				}
				else {
					/* Non-existent element so delta should be non-negative */
					assert(delta >= 0);

					/* Insert sorted */
					if (delta > 0) {
						vect.insert(it, std::move(me));
					}
				}

			}
		};


	};

};