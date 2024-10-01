#include "pch.h"

#include <tuple>
#include <unordered_map>
#include <iostream>

#include "../hgraph/variadic.h"
using namespace hgraph;

TEST(vrdc, hash_tuple) {

	double x = 1.234;
	int64_t n = 1;

	/* Check singletons */
	ASSERT_EQ(std::hash<double>()(x), vrdc::hash_tuple_var<double>()({ x }));
	ASSERT_EQ(std::hash<double*>()(&x), vrdc::hash_tuple_var<double*>()({ &x }));
	ASSERT_EQ(std::hash<int64_t>()(n), vrdc::hash_tuple_var<int64_t>()({ n }));
	ASSERT_EQ(std::hash<int64_t*>()(&n), vrdc::hash_tuple_var<int64_t*>()({ &n }));

	/* Combinations */
	size_t h;
	h = vrdc::hash_tuple_var< double, int64_t >()({ x, n });
	h = vrdc::hash_tuple_var< double*, int64_t* >()({ &x, &n });
	h = vrdc::hash_tuple_var< double*, int64_t*, double >()({ &x, &n, x });

	/* Tuple argument version */
	using tup = std::tuple<int, char, double>;
	tup key = { 2, 'z', 1.23 };
	ASSERT_EQ(
		(vrdc::hash_tuple_var< int, char, double >()(key)),
		(vrdc::hash_tuple< tup >()(key))
	);

	/* Use in map */
	std::unordered_map< tup, size_t, vrdc::hash_tuple<tup> > m;
	size_t val = 1234;
	m.insert({ key, val });
	ASSERT_EQ(val, m.find(key)->second);
}

TEST(vrdc, has_type) {
	using tup = std::tuple<int, char, double>;
	static_assert(vrdc::has_type<int, tup>::value);
	static_assert(!vrdc::has_type<size_t, tup>::value);
}

TEST(vrdc, index_of) {
	static_assert(vrdc::index_of<double, int, char, double>::value == 2);
}

TEST(vrdc, foreach) {
	std::tuple<int, int, int> t{ 1,2,3 };
	int sum;

	/* Lambda capture */
	sum = 0;
	vrdc::foreach(t, [&](int x) { sum += x; });
	ASSERT_EQ(sum, 6);

	/* Argument passing */
	sum = 0;
	vrdc::foreach(t, [](int x, int& sum) { sum += x; }, sum);
	ASSERT_EQ(sum, 6);
}

TEST(vrdc, foreach_if) {
	std::tuple<int, char, int> u{ 1, 'z', 3 };
	int sum;

	sum = 0;
	vrdc::foreach_if<std::is_same, int>(u, [&sum](int x) { sum += x; });
	ASSERT_EQ(sum, 4);

	char z = 0;
	vrdc::foreach_if<std::is_same, char>(u, [&z](char x) { z = x; });
	ASSERT_EQ(z, 'z');
}


template<typename... Ts>
struct Stuff { };


TEST(vrdc, repack_as) {
	using MyTuple = std::tuple<int, char, double>;

	/* Repack as_is */
	using MyTest = vrdc::repack_as<Stuff, MyTuple>::type;
	static_assert(std::is_same_v < MyTest, Stuff<int, char, double> >);

	using MyTestShorthand = vrdc::repack_as_t<Stuff, MyTuple>;
	static_assert(std::is_same_v < MyTestShorthand, Stuff<int, char, double> >);
}

TEST(vrdc, append) {
	using MyTuple = std::tuple<int, char, double>;

	/* Check generic append */
	using MyBiggerVTC = typename vrdc::append_type< double, std::tuple<int, char> >::type;
	static_assert(std::is_same_v<MyTuple, MyBiggerVTC>);
}


/* Predicate: type is double */
template<typename T>
struct IsDouble : std::is_same<T, double> {};

TEST(vrdc, repack_if) {
	/* Repack variadic arguments that satisfy a predicate */
	using MyTuple = std::tuple<int, char, double>;

	/* Check predicate */
	static_assert(IsDouble<double>::value);
	static_assert(!IsDouble<int>::value);

	/* Check conditional tuple repacking */
	using MyTest = typename vrdc::repack_if<IsDouble, MyTuple>::type;
	static_assert(std::is_same_v< MyTest, std::tuple<double> >);
}

/* repack_if recurses too deeply, leads to C1060 error (out of heap space */

/* new implementation, tuple-only */
// see https://stackoverflow.com/questions/44935075/sequence-of-indices-of-tuple-elements-satifying-predicate-in-hana
// see https://codereview.stackexchange.com/questions/115740/filtering-variadic-template-arguments




TEST(vrdc, repack_if2) {
	/* Check using empty tuple */
	using empty_ints = typename vrdc::indices_if<IsDouble, std::tuple<>>::type;

	/* Repack variadic arguments that satisfy a predicate */
	using MyTuple = std::tuple<int, double, char, double>;

	/* Check predicate */
	static_assert(IsDouble<double>::value);
	static_assert(!IsDouble<int>::value);

	using ints = typename vrdc::indices_if<IsDouble, MyTuple>::type;
	static_assert(std::is_same_v< ints, std::index_sequence<1, 3> >);

	using MySubTuple = typename vrdc::subtuple<MyTuple, ints>::type;
	static_assert(std::is_same_v < MySubTuple, std::tuple<double, double> >);

	using MySubTuple2 = typename vrdc::repack_if<IsDouble, MyTuple>::type;
	static_assert(std::is_same_v < MySubTuple2, std::tuple<double, double> >);


	using MyBigTuple = std::tuple<
		//int, char, double,
		//int, char, double,
		//int, char, double,
		//int, char, double,
		//int, char, double,
		//int, char, double,
		//int, char, double,
		//int, char, double,
		//int, char, double,
		//int, char, double,
		//int, char, double,
		//int, char, double,
		//int, char, double,
		//int, char, double,
		//int, char, double,
		//int, char, double,
		//int, char, double,
		//int, char, double,
		//int, char, double,
		//int, char, double,
		//int, char, double,
		//int, char, double,
		//int, char, double,
		//int, char, double,
		//int, char, double,
		//int, char, double,
		//int, char, double,
		//int, char, double,
		//int, char, double,
		//int, char, double,
		//int, char, double,
		//int, char, double,
		//int, char, double,
		//int, char, double,
		//int, char, double,
		//int, char, double,
		//int, char, double,
		//int, char, double,
		int, char, double,
		int, char, double,
		int, char, double,
		int, char, double
	> ;
	using ints2 = typename vrdc::indices_if<IsDouble, MyBigTuple>::type;
	//static_assert(std::is_same_v< ints2, std::index_sequence<2, 5, 8, 11> >);
	//using MySubTuple3 = typename vrdc::repack_if2<IsDouble, MyBigTuple>::type;
}

