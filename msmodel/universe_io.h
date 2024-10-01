#pragma once
#include <string>
#include <iostream>
#include "../msgraph/msgraph_io.h"
#include "universe.h"

/* Defines meta-information about model classes */
struct UniverseReflect {

	/* Construct header for fixed-size array by repeating name with index as suffix */
	static string suffix(const string x, const size_t n) {
		string res;
		res += x + "." + "0";
		for (size_t i = 1; i < n; i++) {
			res += "\t" + x + "." + to_string(i);
		}
		return res;
	}

	/* Graph reference type */
	using graph_type = Universe;

	/* Base struct for trivial (empty) data types. Note that get_trivial::name is undefined to catch errors early */
	struct get_trivial {
		static auto header_tuple() {
			return std::make_tuple();
		}
		template<typename _TrivialData>
		static auto as_tuple(_TrivialData& x) {
			return std::make_tuple();
		}
		using input_tuple_t = std::tuple<>;
	};

	/* Generic access to serialization info for data types. Default (undefined type) prevents serialization. */
	template< typename T >
	struct get {
		using input_tuple_t = std::false_type;
	};

	/* Declarations */

	/****************************************************************************/
	/* MODEL STATE  */
	struct ModelState {
		size_t step;
		ModelState(const size_t step = 0) : step{ step } {}
	};

	template<> struct get< ModelState > {
		static constexpr const char name[] = "ModelState";
		static auto header_tuple() {
			return std::make_tuple("step");
		}
		static auto as_tuple(const ModelState& x) {
			return std::make_tuple(x.step);
		}
		using input_tuple_t = std::tuple<>;
	};

	/****************************************************************************/
	/* GENE */
	template<> struct get< GeneArchetype > {
		static constexpr const char name[] = "GeneArchetype";
		static auto header_tuple() {
			return std::make_tuple(suffix("susceptibility", pressure_base::n_stress_types), "fitness");
		}
		static auto as_tuple(GeneArchetype& x) {
			return std::make_tuple(x.susceptibility, x.fitness);
		}
		using input_tuple_t = std::tuple<array<double, pressure_base::n_stress_types>, double>;
	};

	template<> struct get< Gene > : get_trivial {
		static constexpr const char name[] = "Gene";
	};

	/****************************************************************************/
	/* TRANSPOSON */
	//template<> struct get< TransposonArchetype > {
	//	static constexpr const char name[] = "TransposonArchetype";
	//	static auto header_tuple() {
	//		return std::make_tuple("move");
	//	}
	//	static auto as_tuple(TransposonArchetype& x) {
	//		return std::make_tuple(x.move);
	//	}
	//	using input_tuple_t = std::tuple<double>;
	//};

	//template<> struct get< Transposon > : get_trivial {
	//	static constexpr const char name[] = "Transposon";
	//};

	/****************************************************************************/
	/* PLASMID */
	template<> struct get< PlasmidArchetype > {
		static constexpr const char name[] = "PlasmidArchetype";
		static auto header_tuple() {
			return std::make_tuple("loss", "transfer", "max_count", "fitness");
		}
		static auto as_tuple(PlasmidArchetype& x) {
			return std::make_tuple(x.loss, x.transfer, x.max_count, x.fitness);
		}
		using input_tuple_t = std::tuple<double, double, size_t, double>;
	};

	template<> struct get< Plasmid > : get_trivial {
		static constexpr const char name[] = "Plasmid";
	};

	/****************************************************************************/
	/* CHROMOSOME */

	template<> struct get< ChromosomeArchetype > {
		static constexpr const char name[] = "ChromosomeArchetype";
		static auto header_tuple() {
			return std::make_tuple("fitness", "survival");
		}
		static auto as_tuple(ChromosomeArchetype& x) {
			return std::make_tuple(x.fitness, x.survival);
		}
		using input_tuple_t = std::tuple<double, double>;
	};

	template<> struct get< Chromosome > : get_trivial {
		static constexpr const char name[] = "Chromosome";
	};

	/****************************************************************************/
	/* CELL */
	template<> struct get< CellArchetype > {
		static constexpr const char name[] = "CellArchetype";
		static auto header_tuple() {
			return std::make_tuple();
		}
		static auto as_tuple(CellArchetype& x) {
			return std::make_tuple();
		}
		using input_tuple_t = std::tuple<>;
	};

	template<> struct get< Cell > : get_trivial {
		static constexpr const char name[] = "Cell";
	};

	/****************************************************************************/
	/* PATCH */
	template<> struct get< PatchArchetype > {
		static constexpr const char name[] = "PatchArchetype";
		static auto header_tuple() {
			return std::make_tuple("capacity", suffix("pressure", pressure_base::n_stress_types));
		}
		static auto as_tuple(PatchArchetype& x) {
			return std::make_tuple(x.capacity, x.pressure);
		}
		using input_tuple_t = std::tuple<double, array<double, pressure_base::n_stress_types>>;
	};

	template<> struct get< Patch > : get_trivial {
		static constexpr const char name[] = "Patch";
	};

	/****************************************************************************/
	/* EMPTY EDGE */
	template<> struct get< std::monostate > : get_trivial {
		static constexpr const char name[] = "NoData";
	};

	/****************************************************************************/
	/* MULTIPLICITY */
	template<> struct get< Multiplicity > {
		static constexpr const char name[] = "Multiplicity";
		static auto header_tuple() {
			return std::make_tuple("multiplicity");
		}
		static auto as_tuple(Multiplicity& x) {
			return std::make_tuple(x.mult);
		}
		using input_tuple_t = std::tuple< size_t >;
	};

	/****************************************************************************/
	/* PLASMID COMPATIBILITY RANGE */

	template<> struct get< PlasmidRange > : get_trivial {
		static constexpr const char name[] = "PlasmidRange";
	};

	/****************************************************************************/
	/* DIFFUSION */
	template<> struct get< Diffusion > {
		static constexpr const char name[] = "Diffusion";
		static auto header_tuple() {
			return std::make_tuple("diffusion");
		}
		static auto as_tuple(Diffusion& x) {
			return std::make_tuple(x.diff);
		}
		using input_tuple_t = std::tuple< double >;
	};
};

struct UniverseIO : hgraph_io< UniverseReflect > {};