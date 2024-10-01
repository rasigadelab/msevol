#pragma once

#include "msgraph_decl.h"
#include "meta.h"
#include "multiset.h"

namespace hgraph {

	template<typename _Cfg>
	class msgraph : public graph< _Cfg >
	{
	public:
		/**********************************************/
		/* TYPE ACCESSORS */
		
		/* Declared vertex types */
		using vertex_types = typename _Cfg::vertex_types;
		/* Declared edge types */
		using edge_types = typename _Cfg::edge_types;
		/* Base graph */
		using base_graph = graph<_Cfg>;
		/* Vertex */
		template<typename _VData>
		using V = typename base_graph::template V< _VData >;
		/* Edge */
		template<typename _Source, typename _Target, typename _EData>
		using E = typename base_graph::template E< _Source, _Target, _EData >;

	public:

		/* Remove vertex from graph and forget alterations/barcodes if applicable. */
		template< typename _VData >
		void erase_vertex(V< _VData >& element) {
			if constexpr (std::is_base_of_v<population_container, _VData>) {
				forget_alteration(element);
				//forget_barcode(element);
				mset_forget(element);
			}
			base_graph::erase_vertex(element);
		}

		/**********************************************/
		/* MULTIPLICITY MANAGEMENT */

		/* Propagate count changes to all sub-elements */
		template<typename _VData>
		void shift_count(V< _VData >& x, const int64_t delta) {
			static_assert(std::is_base_of_v<containee, _VData>, "x must be a containee");

			/* Can't remove more than multiplicity (just a safeguard)*/
			assert(x.data.count() + delta >= 0);
			// FIXME Additional temp safeguards
			//assert(x.data.count() < 1e7);
			//assert(delta < 1e7);
			//assert(delta > -1e7);

			/* Update local count*/
			x.data.count_ += delta;

			/* Lambda tail recursion: for each incoming edge list, compute count change
				as the current count change times multiplicity */
			auto propagateCounts = [this](auto& list, const int64_t delta_current) {
				for (auto& edge : list) {
					const int64_t next_delta = delta_current * edge.data.mult;
					if (next_delta != 0) shift_count(*edge.source, next_delta);
				}
			};

			this->template foreach_edge_type<Mlt>(
				x.edges_in_tuple,
				propagateCounts, delta
				);

			/* Delete vertex ? */
			if constexpr (std::is_base_of_v<container, _VData>) {
				if (delta < 0 && x.data.count() == 0 && x.active()) {
					erase_vertex(x);
				}
			}
		}

		/* SHOULD RETURN VOID Update multiplicity and propagate to upper scale. Argument can be negative */
		template<typename _Source, typename _Target>
		E< _Source, _Target, Mlt >&
			shift_mult(E< _Source, _Target, Mlt >& e, const int64_t delta) {

			// FIXME temp safeguards
			assert(delta < 1e7);
			assert(delta > -1e7);

			/* Propagation is limited to immediate upper scale only, hence the
			target element should not be assigned */
			/* Using target.count fails because of Patch static count */
			/* Need to check true outdegree */
			/* WIP: Patch is a singleton, might alleviate the problem using static check ? */
			size_t outdegree = 0;
			this->template foreach_edge_type< Mlt >(e.target->edges_out_tuple,
				[&outdegree](auto& list) {
					//for (auto edge = m.first(); edge != nullptr; m.increment(edge)) {
					for (auto& edge : list) {
						outdegree += edge.data.mult;
					}
				});
			assert((outdegree == 0) && "cannot update multiplicity of already assigned element");

			/* Update edge data */
			assert((int64_t)e.data.mult + delta >= 0);
			e.data.mult += delta;
			assert(e.data.mult >= 0);

			/* Update count */
			const int64_t delta_count = delta * e.target->data.count();
			if (delta_count != 0) shift_count(*e.source, delta_count);

			/* Launch appropriate callback */			
			_Cfg::template callback<_Source, _Target>()(e.source->data, e.target->data, delta);

			/* Delete edge ? */
			if (delta < 0 && e.data.mult == 0 && e.active()) {
				this->erase_edge(e);
			}
			return e;
		};
		
		/**********************************************/
		/* VERTEX ASSIGNATION */

		/* Assign containees to a container, return edge reference. */
		template<typename _Source, typename _Target>
		E< _Source, _Target, Mlt >&
		assign(V< _Source >& source, V< _Target >& target, const size_t mult = 1) {

			assert(mult != 0);

			/* No self-edge for assignation */
			if constexpr (std::is_same_v< _Source, _Target >) {
				assert(&source != &target);
			}

			E< _Source, _Target, Mlt >* e = this->template edge_ptr<Mlt>(source, target);
			if (e == nullptr) e = this->insert_edge_ptr(source, target, Mlt{});

			shift_mult(*e, mult);
			return *e;
		}

		/* Assign containees to a container, return edge reference. */
		template<typename _Source, typename _Target>
		E< _Source, _Target, Mlt >&
		assign(const size_t source_id, const size_t target_id, const size_t mult = 1) {

			assert(mult != 0);
			V< _Source >& source = this->template vertex< _Source >(source_id);
			V< _Target >& target = this->template vertex< _Target >(target_id);

			return assign(source, target, mult);
		}

		/* Attach an element to another with optional multiplicity */
		template<typename _Source, typename _Target>
		std::shared_ptr< E< _Source, _Target, Mlt > >
		assign_shared_ptr(V< _Source >& source, V< _Target >& target, const size_t mult = 1) {

			assert(mult != 0);

			/* No self-edge for assignation */
			if constexpr (std::is_same_v< _Source, _Target >) {
				assert(&source != &target);
			}

			auto e = this->template edge_shared_ptr<Mlt>(source, target);
			if (e == nullptr) e = this->insert_edge_shared_ptr(source, target, Mlt{});

			shift_mult(*e, mult);
			return e;
		}

		template<typename _Source, typename _Target>
		std::shared_ptr< E< _Source, _Target, Mlt > >
		assign_shared_ptr(const size_t source_id, const size_t target_id, const size_t mult = 1) {
			V< _Source >& source = this->template vertex< _Source >(source_id);
			V< _Target >& target = this->template vertex< _Target >(target_id);
			return assign_shared_ptr(source, target, mult);
		}

		/******************************************************************/
		/* ALTERATION KEY */

		/* Edge-like structure that serves as hashmap key for alterations */
		template<typename _Source, typename _Target>
		struct alteration_key {
			static_assert(std::is_base_of_v<containee, _Source>, "Source must be a containee");
			static_assert(std::is_base_of_v<container, _Target>, "Target must be a container");
			V< _Source >* const source;
			V< _Target >* const target;
			const int64_t delta;

			explicit alteration_key(V< _Source >& source_, V< _Target >& target_, const int64_t delta_) :
				source{ &source_ }, target{ &target_ }, delta{ delta_ }
			{
				assert(delta != 0);
			}

			struct hash {
				size_t operator()(const alteration_key<_Source, _Target>& v) const {
					return	(std::hash<const V< _Source >*>()(v.source) << 1) ^
						(std::hash<const V< _Target >*>()(v.target) << 3) ^
						(std::hash<int64_t>()(v.delta) << 5);
				}
			};

			inline bool operator==(const alteration_key<_Source, _Target>& x) const {
				return	source == x.source &&
					target == x.target &&
					delta == x.delta;
			}

			inline bool operator!=(const alteration_key<_Source, _Target>& x) const {
				return !(x == *this);
			}
		};

		/******************************************************************/
		/* ALTERATION MAP */

	private:
		/* Map of alterations */
		template<typename _Source, typename _Target>
		using as_alteration_map = std::unordered_map<
			const alteration_key<_Source, _Target>, V< _Target >*,
			typename alteration_key<_Source, _Target>::hash
		>;
		/* Tuple of alteration maps */
		using alteration_map_tuple_t = typename meta<base_graph>::
			template transform_edges< as_alteration_map, typename meta<base_graph>::multiplicity_edge_types>::type;
		alteration_map_tuple_t alteration_map_tuple;

	public:
		/* Access edge type-specific alteration map */
		template<typename _Source, typename _Target>
		inline as_alteration_map<_Source, _Target>& alteration_map() {
			return std::get< as_alteration_map<_Source, _Target> >(alteration_map_tuple);
		}

		/******************************************************************/
		/* VARIANT MULTIMAP */

	private:
		/* Map of alteration's resulting variant, needed to remove alterations to a dead variant */
		template<typename _Source, typename _Target>
		using as_variant_map = std::unordered_multimap<
			const V< _Target >*, const alteration_key<_Source, _Target>
		>;
		/* Tuple of variant multimaps */
		using variant_map_tuple_t = typename meta<base_graph>::
			template transform_edges< as_variant_map, typename meta<base_graph>::multiplicity_edge_types>::type;
		variant_map_tuple_t variant_map_tuple;

	public:
		/* Access edge type-specific variant multimap */
		template<typename _Source, typename _Target>
		inline as_variant_map<_Source, _Target>& variant_map() {
			return std::get< as_variant_map<_Source, _Target> >(variant_map_tuple);
		}

		/******************************************************************/
		/* ALTERATION LIFECYCLE */

				/* Register an alteration */
		template<typename _Source, typename _Target>
		inline bool register_alteration(V< _Source >& s, V< _Target >& from, V< _Target >& to, const int64_t delta) {
			static_assert(std::is_base_of_v<container, _Target>);

			auto& altmap = alteration_map<_Source, _Target>();
			auto& varmap = variant_map<_Source, _Target>();

			/* Register forward alteration */
			alteration_key< _Source, _Target > key_f(s, from, delta);
			altmap.insert({ key_f, &to });
			varmap.insert({ &to, key_f });
			varmap.insert({ &from, key_f });

			/* Register backward alteration */
			alteration_key< _Source, _Target > key_r(s, to, -delta);
			altmap.insert({ key_r, &from });
			varmap.insert({ &to, key_r });
			varmap.insert({ &from, key_r });

			return true;
		}

		/* Kill a variant, associated alterations and barcode */
		template<size_t I = 0, typename _VData>
		inline void forget_alteration(V< _VData >& x) {
			static_assert(std::is_base_of_v<container, _VData>);
			/* Select incoming multiplicity edges */
			using mult_edges_t = typename vrdc::repack_if<
				typename meta<base_graph>::is_multiplicity_edge,
				typename V< _VData >::edges_in_tuple_t>::type;

			if constexpr (std::tuple_size_v< mult_edges_t > > 0) {
				/* Select source type of Ith edge type */
				using source_t = typename std::tuple_element_t<I, mult_edges_t>::source_t;
				/* Remove from alteration and variant maps */
				auto& varmap = variant_map<source_t, _VData>();
				auto range = varmap.equal_range(&x);
				for (auto it = range.first; it != range.second;) {
					/* Remove alteration pointing to the variant */
					alteration_map< source_t, _VData >().erase(it->second);
					/* Remove variant's entry pointing to the alteration */
					varmap.erase(it++);
				}
				/* Recurse if required */
				if constexpr ((I + 1) < std::tuple_size_v< mult_edges_t >) {
					return forget_alteration< I + 1 >(x);
				}
			}
		}

		/****************************************************************************************/
		/* Generic multiset barcoding + event dispatching routines */

		public:

			/* Multiset shorthands */
			template<typename _VData>
			using mset_vector = typename multiset< base_graph >::template mset_vector< _VData>;
			template<typename _VData>
			using mset = typename multiset< base_graph >::template mset< _VData>;

			/******************************************************************/
			/* MULTISET/BARCODE MAP */

			/* Predicate: vertex is a population container so it can have a multiset barcode. */
			template<typename _VData> struct is_population_container :
				std::is_base_of<population_container, _VData> {};

			/* Tuple of such population container types */
			using population_container_vertices = typename vrdc::repack_if< is_population_container, vertex_types >::type;

			template<typename _VData>
			struct as_mset_map {
				using type = std::unordered_map<mset< _VData >, V< _VData >*, typename mset< _VData >::hash>;
			};

			using mset_map_tuple_t = typename vrdc::transform< as_mset_map, population_container_vertices>::type;
			mset_map_tuple_t mset_map_tuple;

			/* Access to multiset map for a given vertex type */
			template<typename _VData>
			inline typename as_mset_map<_VData>::type& mset_map() {
				static_assert(std::is_base_of_v< population_container, _VData>);
				return std::get< typename as_mset_map<_VData>::type >(mset_map_tuple);
			}

			/******************************************************************/
			/* MULTISET LIFECYCLE */

			/* Compute multiset and register in multiset map. Return true if successful or false
			if multiset was already present in map. */
			template<typename _VData>
			inline bool mset_register(V< _VData >& x) {
				static_assert(std::is_base_of_v< population_container, _VData>);
				/* Compute the barcode */
				auto ms = mset< _VData >(x);
				/* Insert in barcode map */
				auto ret = mset_map<_VData>().insert({ ms, &x });
				return ret.second;
			}

			/* Register multisets of all existing vertices (typically upon startup) */
			template<size_t I = 0>
			inline void mset_register_all() {
				if constexpr (I < std::tuple_size_v< population_container_vertices >) {
					using data_t = std::tuple_element_t<I, population_container_vertices>;
					for (auto& vertex : this->template vertices< data_t >()) {
						/* Don't check uniqueness here for now */
						mset_register(vertex);
					}
					mset_register_all< I + 1 >();
				}
			}
			/* Remove multiset from map, typically before destructing to container vertex. */
			template<typename _VData>
			inline void mset_forget(V< _VData >& x) {
				/* Remove from barcode map */
				size_t ret = mset_map<_VData>().erase(mset< _VData >(x));
				/* FIXME don't check for barcode existence as long as
				container uniqueness is not enforced at startup */
			}
			/* Create a variant based on multiset barcode */
			template<typename _VData>
			V< _VData >* mset_install_vertex_ptr(mset< _VData >& bc) {

				auto vertex_ptr = this->template insert_vertex_ptr< _VData >();

				/* Scan mset and install containees */
				vrdc::foreach(bc.mset_vector_tuple,
					[this, vertex_ptr](auto& vect) {
						for (auto& it : vect) {
							assign(*it.vertex_ptr, *vertex_ptr, it.mult);
						}
					}
				);

				/* Insert in barcode map */
				auto ret = mset_map<_VData>().insert({ bc, vertex_ptr });

				return vertex_ptr;
			}
			/* Find or create a variant based on a multiset barcode. */
			template<typename _VData>
			V< _VData >* mset_get_vertex_ptr(mset< _VData >& bc) {

				/* Try find vertex in map */
				auto& map = mset_map< _VData >();
				auto it = map.find(bc);
				if (it == map.end()) {
					/* Create vertex */
					return mset_install_vertex_ptr(bc);
				}
				else {
					return it->second;
				}
			}
			/* Find or create a variant based on a multiset barcode. */
			template<typename _VData>
			inline V< _VData >& mset_get_vertex(mset< _VData >& bc) {
				return *mset_get_vertex_ptr(bc);
			}
	}; // class msgraph
}