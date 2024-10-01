#include "pch.h"

#include "../hgraph/typelist.h"
#include <vector>

using namespace hgraph;


/* Metafunction: change a type to a vector */
template<typename T>
struct as_vector {
	using type = std::vector< T >;
};

/* Predicate */
template<typename T> struct IsChar : std::is_same<T, char> {};

TEST(typelist, test) {
	using List = tl<int, double, char>;

	using List2 = typelist<tl, List>;
	static_assert(List2::size == 3);
	static_assert(std::is_same_v< typename List2::at<0>::type, int >);
	static_assert(std::is_same_v< typename List2::at<1>::type, double >);
	static_assert(std::is_same_v< typename List2::at<2>::type, char >);
	static_assert(std::is_same_v< typename List2::at<3>::type, void >);

	using List3 = typelist<std::tuple, std::tuple<int, double, char> >;

	using ListVec = typename typelist<tl, List>::transform< as_vector >::type;
	static_assert(std::is_same_v< ListVec, tl< std::vector<int>, std::vector<double>, std::vector<char> > >);

	using ListTuple = typename typelist<tl, List>::repack_as<std::tuple>::type;
	static_assert(std::is_same_v< ListTuple, std::tuple<int, double, char> >);

	/* Try chain operator */
	using List4 = typename List2::transform< as_vector >::chain::repack_as< std::tuple >::type;

	/* push_back */
	using List5 = typename List2::push_back< float >::chain;
	static_assert(List5::size == 4);
	static_assert(std::is_same_v< typename List5::at<3>::type, float >);

	/* push front */
	using List6 = typename List2::push_front< float >::chain;
	static_assert(List6::size == 4);
	static_assert(std::is_same_v< typename List6::at<0>::type, float >);

	/* extract */
	using List7 = typename List2::extract< std::index_sequence< 0, 2 > >::chain;
	static_assert(List7::size == 2);
	static_assert(std::is_same_v< typename List7::at<0>::type, int >);
	static_assert(std::is_same_v< typename List7::at<1>::type, char >);

	/* indices satisfying a predicate */
	using List8 = typename List2::indices_of< IsChar >::type;
	static_assert(std::is_same_v< List8, std::index_sequence< 2 > >);

	/* types satisfying a predicate */
	using List9 = typename List2::filter< IsChar >::type;
	static_assert(std::is_same_v< List9, tl< char > >);

	/* check that list contains a type */
	static_assert(List2::has_type<int>::value);
	static_assert(!List2::has_type< size_t >::value);

	/* Shorter wrapper for tuples */
	using List10 = typetuple< std::tuple< int, double, char > >;
	static_assert(List10::size == 3);
}

//template< size_t I = 0, typename... Ts >
//inline void add_one(std::tuple< Ts... >& x) {
//	static constexpr size_t N = sizeof...(Ts);
//	if constexpr (N > 0) {
//		if constexpr (I)
//			std::get< I >(x) += 1;
//	}
//
//
//}

/* Metafunction as a function-returning struct */
template<typename T> struct add_one {
	void operator()(T& x) { x++;  }
};

/* Metafunction as an accumulator */
template<typename T> struct metasum {
	size_t operator()(T& x, size_t accum = 0) { return accum + x; }
};

/* Function templates are not OK as template arguments */
//template<typename T> void add_one2(T& x) { x++; }

/* Build function from metafunction and tuple */
//template< template<typename> typename _Metafunc, typename... Ts > struct meta;

template< template<typename> typename _Metafunc, size_t I = 0, typename... Ts, typename... _Args >
void meta( std::tuple< Ts... >& x, _Args... args ) {
	constexpr size_t N = sizeof...(Ts);
	using T_I = std::tuple_element_t< I, std::tuple< Ts... > >;
	_Metafunc< T_I >()( std::get< I >( x ) );
	if constexpr (I + 1 < N) {
		return meta< _Metafunc, I + 1 >(x);
	}
};


// See https://stackoverflow.com/questions/4697180/template-function-as-a-template-argument

TEST(Metafunc, test) {
	/* Test metafunction application on tuples to preserve type safety in lambda-like functions */
	//using List = std::tuple<int, double, float>;
	//List l = {0, 0, 0};

	std::tuple<int, double, float> l = { };

	/* Apply a metafunction +1 */
	meta< add_one >(l);

	EXPECT_EQ(std::get< 0 >(l), 1);
}
