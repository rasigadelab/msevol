#pragma once

#include "model_vertices.h"
#include "model_edges.h"
#include "../msgraph/msgraph.h"

using namespace hgraph;

/*****************************************************************************************/
/* PROPERTY PROPAGATION CALLBACK */

/* Propagation by assignation */
template<typename T>
inline void propagate_as_constant(T& source_prop, T& target_prop, const int64_t delta) {
	target_prop = source_prop;
}

/* Propagation by assignation */
template<typename T, size_t N>
inline void propagate_array_as_constant(
	std::array<T, N>& source_prop,
	std::array<T, N>& target_prop,
	const int64_t delta
) {
	for (size_t i = 0; i < N; ++i) {
		target_prop[i] = source_prop[i];
	}
}

/* Default propagation method with multiplication/division */
inline void propagate_as_product(double& source_prop, double& target_prop, const int64_t delta) {
	if (delta > 0) {
		target_prop *= std::pow(source_prop, delta);
	}
	else if (delta < 0) {
		target_prop /= std::pow(source_prop, -delta);
	}
}

/* Default propagation method with multiplication/division, for arrays */
template<size_t N>
inline void propagate_array_as_product(
	std::array<double, N>& source_prop,
	std::array<double, N>& target_prop,
	const int64_t delta
) {
	if (delta > 0) {
		for (size_t i = 0; i < N; ++i) {
			target_prop[i] *= std::pow(source_prop[i], delta);
		}
	}
	else if (delta < 0) {
		for (size_t i = 0; i < N; ++i) {
			target_prop[i] /= std::pow(source_prop[i], -delta);
		}
	}
}

/* Property propagation callback */
template<typename _Source, typename _Target>
struct propagate {
	inline void operator()(_Source& source, _Target& target, const int64_t delta) {
		//static_assert(false, "Default propagation functor fallback.");
	}
};

/* ARCHETYPE PROPAGATION */

template<> struct propagate< GeneArchetype, Gene > {
	inline void operator()(GeneArchetype& source, Gene& target, const int64_t delta) {
		/* Fitness & susceptibility */
		assert(delta != 0);
		propagate_as_constant(source.fitness, target.fitness, delta);
		propagate_array_as_constant(source.susceptibility, target.susceptibility, delta);
	}
};

//template<> struct propagate< TransposonArchetype, Transposon > {
//	inline void operator()(TransposonArchetype& source, Transposon& target, const int64_t delta) {
//		/* Movement rate */
//		assert(delta != 0);
//		propagate_as_constant(source.move, target.move, delta);
//	}
//};

template<> struct propagate< ChromosomeArchetype, Chromosome > {
	inline void operator()(ChromosomeArchetype& source, Chromosome& target, const int64_t delta) {
		/* Fitness & susceptibility */
		assert(delta != 0);
		propagate_as_constant(source.fitness, target.fitness, delta);
		propagate_as_constant(source.survival, target.survival, delta);
	}
};

template<> struct propagate< PlasmidArchetype, Plasmid > {
	inline void operator()(PlasmidArchetype& source, Plasmid& target, const int64_t delta) {
		/* Fitness & susceptibility */
		assert(delta != 0);
		propagate_as_constant(source.fitness, target.fitness, delta);
		propagate_as_constant(source.loss, target.loss, delta);
		propagate_as_constant(source.transfer, target.transfer, delta);
		propagate_as_constant(source.max_count, target.max_count, delta);
	}
};

template<> struct propagate< PatchArchetype, Patch > {
	inline void operator()(PatchArchetype& source, Patch& target, const int64_t delta) {
		/* Fitness & susceptibility */
		assert(delta != 0);
		propagate_as_constant(source.capacity, target.capacity, delta);
		propagate_array_as_constant(source.pressure, target.pressure, delta);
	}
};

/* PROPERTY PROPAGATION */

template<> struct propagate< Gene, Chromosome > {
	inline void operator()(Gene& source, Chromosome& target, const int64_t delta) {
		/* Fitness & susceptibility */
		assert(delta != 0);
		propagate_as_product(source.fitness, target.fitness, delta);
		propagate_array_as_product(source.susceptibility, target.susceptibility, delta);
	}
};

template<> struct propagate< Gene, Plasmid > {
	inline void operator()(Gene& source, Plasmid& target, const int64_t delta) {
		/* Fitness & susceptibility */
		assert(delta != 0);
		propagate_as_product(source.fitness, target.fitness, delta);
		propagate_array_as_product(source.susceptibility, target.susceptibility, delta);
	}
};

//template<> struct propagate< Gene, Transposon > {
//	inline void operator()(Gene& source, Transposon& target, const int64_t delta) {
//		/* Fitness & susceptibility */
//		assert(delta != 0);
//		propagate_as_product(source.fitness, target.fitness, delta);
//		propagate_array_as_product(source.susceptibility, target.susceptibility, delta);
//	}
//};

//template<> struct propagate< Transposon, Chromosome > {
//	inline void operator()(Transposon& source, Chromosome& target, const int64_t delta) {
//		/* Fitness & susceptibility */
//		assert(delta != 0);
//		propagate_as_product(source.fitness, target.fitness, delta);
//		propagate_array_as_product(source.susceptibility, target.susceptibility, delta);
//	}
//};

//template<> struct propagate< Transposon, Plasmid > {
//	inline void operator()(Transposon& source, Plasmid& target, const int64_t delta) {
//		/* Fitness & susceptibility */
//		assert(delta != 0);
//		propagate_as_product(source.fitness, target.fitness, delta);
//		propagate_array_as_product(source.susceptibility, target.susceptibility, delta);
//	}
//};

template<> struct propagate< Chromosome, Cell > {
	inline void operator()(Chromosome& source, Cell& target, const int64_t delta) {
		/* Fitness, susceptibility and survival */
		assert(delta != 0);
		propagate_as_product(source.fitness, target.fitness, delta);
		propagate_as_product(source.survival, target.survival, delta);
		propagate_array_as_product(source.susceptibility, target.susceptibility, delta);
	}
};

template<> struct propagate< Plasmid, Cell > {
	inline void operator()(Plasmid& source, Cell& target, const int64_t delta) {
		/* Fitness & susceptibility */
		assert(delta != 0);
		propagate_as_product(source.fitness, target.fitness, delta);
		propagate_array_as_product(source.susceptibility, target.susceptibility, delta);
	}
};

template<> struct propagate< Cell, Patch > {
	inline void operator()(Cell& source, Patch& target, const int64_t delta) {
		/* Population count */
		assert(delta != 0);
		target.popsize += delta;
	}
};

/*****************************************************************************************/
/* MODEL GRAPH SCHEMA */

/* Current I/O implementation makes order important for multiplicity edges as
they should be parsed from containee to container (toposort_like) */

/* Graph schema structure */
struct Cfg {
	
	using vertex_types = std::tuple<
		GeneArchetype, Gene,
		//TransposonArchetype, Transposon,
		PlasmidArchetype, Plasmid,
		ChromosomeArchetype, Chromosome,
		CellArchetype, Cell,
		PatchArchetype, Patch
	>;

	using edge_types = std::tuple <
		/* Archetypes */
		Edge< GeneArchetype, Gene, Mlt >,
		//Edge< TransposonArchetype, Transposon, Mlt >,
		Edge< PlasmidArchetype, Plasmid, Mlt >,
		Edge< ChromosomeArchetype, Chromosome, Mlt >,
		Edge< CellArchetype, Cell, Mlt >,
		Edge< PatchArchetype, Patch, Mlt >,
		/* Containers */
		//Edge< Gene, Transposon, Mlt >,
		Edge< Gene, Plasmid, Mlt >,
		Edge< Gene, Chromosome, Mlt >,
		//Edge< Transposon, Transposon, Mlt >,
		//Edge< Transposon, Plasmid, Mlt >,
		//Edge< Transposon, Chromosome, Mlt >,
		Edge< Plasmid, Cell, Mlt >,
		Edge< Chromosome, Cell, Mlt >,
		Edge< Cell, Patch, Mlt >,
		/* Diffusion between patches */
		Edge< Patch, Patch, Diff >,
		/* Plasmid transfer data */
		Edge< PlasmidArchetype, ChromosomeArchetype, PlasmidRange >,
		Edge< Plasmid, Patch, Conjugation >//,
		/* Gene jump data */
		//Edge< Transposon, Cell, GeneJump >
	>;
	
	template<typename _Source, typename _Target>
	using callback = propagate< _Source, _Target >;
};

class Universe : public msgraph< Cfg > {
public:
	/* Clarify vertex typename */
	template<typename _VData>
	using V = typename graph< Cfg >::template V< _VData >;

	/* Clarify edge typename */
	template<typename _Source, typename _Target, typename _EData>
	using E = typename graph< Cfg >::template E< _Source, _Target, _EData >;

	/* Extract element type from archetype */
	template< typename _Arch >
	using elem_t = typename _Arch::elem_t;

	/* Extract archetype type from element type */
	template< typename _VData >
	using arch_t = typename _VData::arch_t;

	/* Register a new container with a type element */
	template<typename _Arch>
	std::shared_ptr< V< elem_t< _Arch > > >
		insertTypedContainerShared(V< _Arch >& type, const size_t id) {
		auto element = this->template insert_vertex_shared_ptr< elem_t< _Arch > >(id);
		assign(type, *element);
		return element;
	}

	/* Register a new element with a type element  */
	template<typename _Arch>
	inline std::shared_ptr< V< elem_t< _Arch > > >
		insertTypedContainerShared(V< _Arch >& type) {
		return insertTypedContainerShared(type, this->template next_vertex_id< elem_t< _Arch > >());
	}

	/* Register a new container with a type element */
	template<typename _Arch>
	V< elem_t< _Arch > >& insertTypedContainer(V< _Arch >& type, const size_t id) {
		auto& element = this->template insert_vertex< elem_t< _Arch > >(id);
		assign(type, element);
		return element;
	}

	/* Register a new element with a type element  */
	template<typename _Arch>
	inline V< elem_t< _Arch > >& insertTypedContainer(V< _Arch >& type) {
		return insertTypedContainer(type, this->template next_vertex_id< elem_t< _Arch > >());
	}

	/* Find archetype of an element */
	template< typename _VData >
	inline std::shared_ptr< V< arch_t< _VData > > >
		getArchetypeShared (V< _VData >& element) {
		static_assert(std::is_base_of_v<container, _VData>, "x must be a container");
		return element.edges_from< arch_t< _VData >, Mlt >().front().source->shared_from_this();
	}

	/* Find archetype of an element */
	template< typename _VData >
	inline V< arch_t< _VData > >& getArchetype(V< _VData >& element) {
		static_assert(std::is_base_of_v<container, _VData>, "x must be a container");
		return *(element.edges_from< arch_t< _VData >, Mlt >().front().source);
	}
};