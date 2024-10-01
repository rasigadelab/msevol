#pragma once
#include "../hgraph/hgraph.h"

using namespace hgraph;

struct A {};
struct B {};

/* Edge data types */
struct X {};
struct Y {};
struct XY : X, Y {};

/* Declare vertices */
using NodeTypes = std::tuple< A, B >;

/* Declare edges */
using EdgeTypes = std::tuple<
	Edge<A, B, X>, Edge<A, B>, Edge<A, B, Y>, Edge<A, B, XY>
>;

//class Graph : public graph< NodeTypes, EdgeTypes > {};

/* V2: configuration struct */
struct Cfg {
	using vertex_types = std::tuple< A, B >;
	using edge_types = std::tuple<
		Edge<A, B, X>, Edge<A, B>, Edge<A, B, Y>, Edge<A, B, XY>
	>;
};

class Graph : public graph< Cfg > {};