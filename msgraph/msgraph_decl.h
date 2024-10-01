#pragma once

/***************************************************************************/
/* MULTISET HIERARCHY GRAPH */

/*
Special type of heterogeneous graph meant to handle hierarchical structures
that involve multisets, similar to P-systems (membrane computing)

MOTIVATION.
	- represent a hierarchy of heterogeneous multisets that can be manipulated
	with a limited set of graph rewriting functions.
	- enables compact P-system computations where populations of identical individuals
	are modelled as a single vertex.

COUNT PROPAGATION.
	- each change of the multiplicity of an element in a multiset is propagated to all its
	containees, allowing to track their count.

PROPERTY PROPAGATION.
	- properties of a multiset are defined as a function of its content. Each change of the
	multiplicity of a containee modifies recursively the properties of the container.
	- the details of the propagation method can vary with properties.
	- the propagation method is passed as a callback template functor. Callback for
	property propagation (see default_callback) is called each time a multiplicity
	edge is updated. The callback's operator() member takes the data types of source
	and target vertices and the change of multiplicity as arguments to update
	the properties of the target.

MULTISET UNIQUENESS ENFORCEMENT.
	- multisets evolves by changing the multiplicity of their contained elements (containees).
	- identical multisets often arise repeatedly during evolution of the system.
	- to avoid redundant calculations on identical multisets, msgraphs enforces their uniqueness:
	if a change of a multiset results into another existing multiset, then the existing multiset is
	retrieved and taken as the result of the change
	- uniqueness enforcement uses a 2-phase lookup based on alteration tracking and barcoding
	- each multiplicity change that alters a multiset is registered in the form of an 'alteration'. If
	the same change occurs, the alteration is retrieved and the result of the alteration is the
	previous result. Alterations are bidirectional. Alteration tracking is in amortized constant time.
	- if no existing alteration is found, a preexisting result is found by searching if another container
	has the same containees and multiplicities (barcoding). Barcoding complexity is linear in the number
	of containees.
	- if both alteration lookup and barcode lookup fail, the resulting multiset is considered a novel
	container. The new contained is added to the graph and its alterations and barcode are registered.

*/

#include "../hgraph/hgraph.h"
 

namespace hgraph {

/**********************************************/
/* MSGRAPH CONTAINER / CONTAINEE TAGS */

	/*
	template inheritance is verbose, but we need it...
	see https://stackoverflow.com/questions/610245/where-and-why-do-i-have-to-put-the-template-and-typename-keywords
	prefix base functions with this->template
	*/

	/* Tag for container, a vertex representing a multiset of other vertices (the containees) */
	struct container {};
	/* Tag for container that behaves as a population and can have multiplicity > 1*/
	struct population_container : container {};
	/* Tag for container that behaves as an individual and always has multiplicity <= 1 */
	struct singleton_container : container {};

	/* A vertex that can be contained by a container */
	class containee {
		template<typename> friend class msgraph;
	protected:
		size_t count_;
	public:
		containee() : count_{ 0 } {}
		inline size_t count() { return count_; }
	};

	/* Multiplicity edge data structure */
	struct Multiplicity {
		size_t mult;
		explicit Multiplicity(const size_t m = 0) : mult{ m } {}
	};
	using Mlt = Multiplicity;

	/* Default callback for property propagation */
	template<typename _Source, typename _Target>
	struct default_callback {
		inline void operator()(_Source& s, _Target& t, const int64_t m) {}
	};

	//template<typename _Vs, typename _Es, template<typename, typename> typename _Cb = default_callback>
	//class msgraph;

	template<typename T>
	static inline std::vector<T>& operator-=(std::vector<T>& lhs, const std::vector<T>& rhs) {
		assert(lhs.size() == rhs.size());
		for (size_t i = 0; i < lhs.size(); ++i) {
			lhs[i] -= rhs[i];
		}
		return lhs;
	}

} // namespace