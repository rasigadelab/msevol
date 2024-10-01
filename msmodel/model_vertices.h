#pragma once
#include <cmath>
#include "../msgraph/msgraph.h"
using namespace hgraph;

/*****************************************************************************************/
/* VERTEX PROPERTIES */

struct pressure_base {
	static constexpr size_t n_stress_types = 2;
};

struct Susceptibility : pressure_base {
	std::array<double, n_stress_types> susceptibility;
	Susceptibility(const double s = 1.) { susceptibility.fill( s ); }
	Susceptibility(const std::array<double, n_stress_types>& s) : susceptibility{ s } {}
};

struct Fitness {
	double fitness;
	Fitness(const double f = 1.) : fitness{ f } {}
};

struct Survival {
	double survival;
	Survival(const double s = 1.) : survival{ s } {}
};

struct PlasmidMovement {
	double loss;
	double transfer;
	size_t max_count;
	PlasmidMovement(const double loss = 0., const double transfer = 0., const size_t max_count = 1) :
		loss{ loss }, transfer{ transfer }, max_count{ max_count } {}
};

struct Capacity {
	double capacity;
	Capacity(const double c = 1000) : capacity{ c } {}
};

struct Pressure : pressure_base {
	std::array<double, n_stress_types> pressure;
	Pressure(const double p = 0) { pressure.fill(p); }
	Pressure(const std::array<double, n_stress_types>& p) : pressure{ p } {}
};

struct Population {
	size_t popsize;
	Population(const size_t p = 0) : popsize{ p } {}
};

//struct TransposonMovement {
//	double move;
//	TransposonMovement(const double move = 0.) : move{ move } {}
//};

/****************************************************************************/
/* Declare vertex types */

/****************************************************************************/
/* GENE */

class Gene;

class GeneArchetype :
	public Susceptibility,
	public Fitness,
	public containee
{
public:
	/* Access to vertex type */
	using elem_t = Gene;

	GeneArchetype() {};
	GeneArchetype(const std::array<double, n_stress_types>& s, const double f) :
		Susceptibility{ s }, Fitness{ f } {};

	/* DEPRECATED set() interface */
	GeneArchetype& set(const double suscept_, const double fitn_) {
		assert(suscept_ > 0.);
		assert(suscept_ <= 1.);
		assert(fitn_ > 0.);
		assert(fitn_ <= 1.);
		susceptibility.fill(suscept_);
		fitness = fitn_;
		return *this;
	}
	GeneArchetype& set(const std::array<double, pressure_base::n_stress_types> suscept_, const double fitn_) {
		for (size_t i = 0; i < pressure_base::n_stress_types; i++) {
			assert(suscept_[i] > 0.);
			assert(suscept_[i] <= 1.);
		}
		assert(fitn_ > 0.);
		assert(fitn_ <= 1.);
		susceptibility = suscept_;
		fitness = fitn_;
		return *this;
	}
};

class Gene :
	public containee,
	public population_container,
	public Fitness,
	public Susceptibility
{
public:
	using arch_t = GeneArchetype;
};

/****************************************************************************/
/* TRANSPOSON */

//class Transposon;

//class TransposonArchetype :
//	public TransposonMovement,
//	public containee
//{
//public:
//	using elem_t = Transposon;
//
//	TransposonArchetype() {}
//	TransposonArchetype(const double move) : TransposonMovement{ move } {}
//};

//class Transposon :
//	public TransposonMovement,
//	public Susceptibility,
//	public Fitness,
//	public population_container,
//	public containee
//{
//public:
//	using arch_t = TransposonArchetype;
//};


/****************************************************************************/
/* CHROMOSOME */

class Chromosome;

class ChromosomeArchetype :
	public Fitness,
	public Survival,
	public containee
{
public:
	/* Access to vertex type */
	using elem_t = Chromosome;

	ChromosomeArchetype() {}
	ChromosomeArchetype(const double fitness, const double survival) :
		Fitness{ fitness }, Survival{ survival } {}

	ChromosomeArchetype& set(const double fitn_, const double surv_) {
		fitness = fitn_;
		survival = surv_;
		return *this;
	}
};

class Chromosome :
	public containee,
	public population_container,
	public Fitness,
	public Susceptibility,
	public Survival
{
public:
	using arch_t = ChromosomeArchetype;
};

/****************************************************************************/
/* PLASMID */
class Plasmid;

class PlasmidArchetype :
	public PlasmidMovement,
	public Fitness,
	public containee
{
public:
	using elem_t = Plasmid;

	PlasmidArchetype() {}
	PlasmidArchetype(const double loss, const double transfer, const size_t max_count, const double fitness) :
		PlasmidMovement{ loss, transfer, max_count }, Fitness{ fitness } {}

	PlasmidArchetype& set(const double loss_, const double transfer_, const double fitn_) {
		loss = loss_;
		transfer = transfer_;
		fitness = fitn_;
		return *this;
	}
};

class Plasmid :
	public containee,
	public population_container,
	public Fitness,
	public Susceptibility,
	public PlasmidMovement
{
public:
	using arch_t = PlasmidArchetype;
};

/****************************************************************************/
/* CELL */
class Cell;

class CellArchetype :
	public containee
{
public:
	using elem_t = Cell;
};

class Cell :
	public containee,
	public population_container,
	public Fitness,
	public Susceptibility,
	public Survival
{
public:
	using arch_t = CellArchetype;
};

/****************************************************************************/
/* PATCH */
class Patch;

class PatchArchetype :
	public Capacity,
	public Pressure,
	public containee
{
public:
	/* Access to vertex type */
	using elem_t = Patch;

	PatchArchetype() {}
	PatchArchetype(const double capacity, const std::array<double, pressure_base::n_stress_types> pressure) :
		Capacity{ capacity }, Pressure{ pressure } {}

	/* Common-interface setter */
	PatchArchetype& set(const double capacity_, const double pressure_) {
		capacity = capacity_;
		pressure.fill(pressure_);
		return *this;
	}
	PatchArchetype& set(const double capacity_,
		const std::array<double, pressure_base::n_stress_types> pressure_V_) {
		capacity = capacity_;
		pressure = pressure_V_;
		return *this;
	}
};

/* Note initialization of count = 1 because patch is singleton */

class Patch :
	public containee,
	public singleton_container,
	public Capacity,
	public Pressure,
	public Population
{
public:
	using arch_t = PatchArchetype;
	Patch() { count_ = 1; }
};





