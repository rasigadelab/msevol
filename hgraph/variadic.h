#pragma once

#include <type_traits>
#include <tuple>
#include <utility>

namespace hgraph {

	/* Static routines to simplify use of tuples */
	class vrdc {

		vrdc() = delete;

/****************************************************************/
/* TUP_APPLY: FUNCTOR-BASED TUPLE VISITOR */
	public:
		/* Functor-based tuple visitor */
		template<template<typename...> typename _Functor, size_t I = 0, typename _Tuple, typename... _Args>
		static inline void tup_apply(_Tuple& t, _Args&... args) {
			if constexpr (I < std::tuple_size_v< _Tuple >) {
				using T = std::tuple_element_t< I, _Tuple >;
				_Functor< T >()(std::get< I >(t), args...);
				tup_apply< _Functor, I + 1, _Tuple, _Args... >(t, args...);
			}
		}

/****************************************************************/
/* HASH ROUTINES FOR TUPLE-BASED KEY */

		/* Recursive hash computation on tuple - variadic version */
	private:
		template<size_t I, typename... T>
		static inline size_t hash_tuple_var_(const std::tuple<T...>& x) {
			constexpr size_t N = sizeof...(T);
			static_assert(N > 0, "Empty tuple");
			using T_I = std::tuple_element_t< I, std::tuple<T...> >;
			if constexpr (I + 1 < N) {
				return std::hash<T_I>()(std::get<I>(x)) ^ hash_tuple_var_<I + 1, T...>(x);
			}
			if constexpr (I + 1 == N) {
				return std::hash<T_I>()(std::get<I>(x));
			}
		};
	public:
		/* Compute hash on tuple. Use vrdc::hash_tuple<T...>()(tup).
		Pass vrdc::hash_tuple<T...>	as hash argument. */
		template<typename... T>
		struct hash_tuple_var {
			size_t operator()(const std::tuple<T...>& x) const
			{
				size_t h = hash_tuple_var_<0, T...>(x);
				return  h;
			}
		};

		/* Recursive hash computation on tuple - tuple argument version */
	private:
		template<size_t I, typename T>
		static inline size_t hash_tuple_(const T& x) {
			constexpr size_t N = std::tuple_size_v<T>;
			static_assert(N > 0, "Empty tuple");
			using T_I = std::tuple_element_t< I, T >;
			if constexpr (I + 1 < N) {
				return std::hash<T_I>()(std::get<I>(x)) ^ hash_tuple_<I + 1, T>(x);
			}
			if constexpr (I + 1 == N) {
				return std::hash<T_I>()(std::get<I>(x));
			}
		};
	public:
		/* Compute hash on tuple. Use vrdc::hash_tuple< tuple<T...> >()(tup).
		Pass vrdc::hash_tuple< tuple<T...> > as hash argument. */
		template<typename T>
		struct hash_tuple {
			size_t operator()(const T& x) const
			{
				size_t h = hash_tuple_<0, T>(x);
				return  h;
			}
		};

/****************************************************************/
/* HAS_TYPE: SEARCH FOR TYPE IN TUPLE */
	public:
		/* Check whether a tuple contains a type */
		template <typename T, typename Tuple>
		struct has_type;
		/* Check whether a tuple contains a type */
		template <typename T, typename... Us>
		struct has_type<T, std::tuple<Us...>> : std::disjunction<std::is_same<T, Us>...> {};

		/****************************************************************/
		/* INDEX_OF: FIND INDEX OF TYPE IN PARAMETER PACK */
	private:
		template<size_t I, typename T, typename ...Ts>
		struct index_of_;

		template<size_t I, typename T>
		struct index_of_<I, T> {
			static constexpr size_t value = I;
		};

		template<size_t I, typename T1, typename T2, typename ...Ts>
		struct index_of_<I, T1, T2, Ts...> {
			static constexpr size_t value =
				std::is_same<T1, T2>::value ?
				I : index_of_<I + 1, T1, Ts...>::value;
		};

	public:
		/* Find index of typename in parameter pack */
		template<typename T1, typename T2, typename ...Ts>
		struct index_of : index_of_<0, T1, T2, Ts...> {};

/****************************************************************/
/* INDEX_IN_TUPLE: Find index of typename in tuple */

			/*Find index of typename in tuple */
		template<typename T, typename... Ts>
		struct index_in_tuple;
		/*Find index of typename in tuple */
		template<typename T, typename... Ts>
		struct index_in_tuple<T, std::tuple<Ts...>> : index_of<T, Ts...> {};

/****************************************************************/
/* FOREACH: Apply function on tuple */

			/* Recursively apply function on tuple elements and pass additional arguments.
			Wrap methods in lambda expression. */
		template<size_t I = 0, typename... T, typename F, typename... A>
		static inline void foreach(std::tuple<T...>& t, F func, A&... args) {
			if constexpr (sizeof...(T) > 0) {
				func(std::get< I >(t), args...);
				if constexpr ((I + 1) < sizeof...(T)) {
					foreach< I + 1 >(t, func, args...);
				}
			}
		}

/****************************************************************/
/* FOREACHIF: Apply function on tuple conditional on predicate */

		/* Calls Predicate< type, Pargs...> on each tuple element type */
	private:
		template< size_t I, template<typename...> typename Predicate,
			typename... Pargs, typename... T, typename F, typename... A>
			static inline void foreach_if_impl(std::tuple<T...>& t, F func, A&... args) {
			if constexpr (sizeof...(T) > 0) {
				using T_I = std::tuple_element_t< I, std::tuple<T...> >;
				if constexpr (
					Predicate< T_I, Pargs... >::value
					) {
					func(std::get< I >(t), args...);
				}
				if constexpr ((I + 1) < sizeof...(T)) {
					foreach_if_impl< I + 1, Predicate, Pargs... >(t, func, args...);
				}
			}
		}

	public:
		/* Apply function on tuple elements whose type satisfy a predicate. */
		template< template<typename...> typename Predicate,
			typename... Pargs, typename... T, typename F, typename... A>
			static inline void foreach_if(std::tuple<T...>& t, F func, A&... args) {
			return foreach_if_impl<0, Predicate, Pargs...>(t, func, args...);
		}

/********************************************************************/
/* VARIADIC TYPE MANIPUATION */

	public:
		/* Variadic type manipulation: transform _From<T...> into _To<T...> */
		template<template<typename...> class T, typename>
		struct repack_as { };

		template<
			template<typename...> typename _To,
			template<typename...> typename _From,
			typename... Ts
		>
			struct repack_as<_To, _From<Ts...>>
		{
			using type = _To<Ts...>;
		};

		/* Variadic type manipulation: transform _From<T...> into _To<T...> */
		template< template<typename...> typename _To, typename _From >
		using repack_as_t = typename repack_as<_To, _From>::type;

	public:
		/* Variadic type manipulation: append a type to a variadic type */
		template<typename, typename>
		struct append_type;

		/* Variadic type manipulation: append a type to a variadic type */
		template<typename T, template<typename...> typename _VTC, typename... Ts>
		struct append_type< T, _VTC<Ts...> > {
			using type = _VTC< Ts..., T >;
		};

		/* Variadic type manipulation: build tuple from integer sequence */
		template<typename, typename>
		struct subtuple;

		/* Variadic type manipulation: build tuple from integer sequence */
		template<typename _Tuple, size_t... Is>
		struct subtuple<_Tuple, std::index_sequence<Is...> > {
			using type = std::tuple<typename std::tuple_element<Is, _Tuple>::type...>;
		};

		/* Variadic type manipulation: apply a transformation to all tuple elements. Transformation
		should not be nested. */
		template<template<typename> typename, typename> struct transform;
		template<template<typename> typename _Func, typename... T>
		struct transform< _Func, std::tuple< T... > > {
			using type = std::tuple< typename _Func< T >::type... >;
		};

		/* Variadic type manipulation: apply a transformation to all tuple elements. Transformation
		should not be nested. */
		template<template<typename> typename, typename> struct transform_as;
		template<template<typename> typename _Func, typename... T>
		struct transform_as< _Func, std::tuple< T... > > {
			using type = std::tuple< _Func< T >... >;
		};

/********************************************************************/
/* INDICES_IF: Indices of tuple types satisfying a predicate */

	public:

		/* Return the sequence of types indices that satisfy a predicate */
	template <template<typename> typename _Pred, typename _Tuple>
	class indices_if {

		/***************************
		* Implement compile-time control sequence on types.
		* For type index I, recurse< I > determines whether the predicate is satisfied
		* and calls increment<true> or increment<false> to add or not the index to the list.
		* In turn, increment< I > calls recurse< I + 1 > to continue iteration until
		* the base case I == N is found.

		* This less-recursive implementation is a workaround for the naive std::conditional_t<> recursion
		* which leads to 2^n exponential compilation time if recursive calls are made
		* in both true and false condition branches, as all branches are evaluated recursively.
		*
		* Note that current implementation is quadratic in the number M of types satisfying
		* the predicate because M types need to be deduced until the sequence is returned.

		* see http://ericniebler.com/2014/11/13/tiny-metaprogramming-library/
		* see https://softwareengineering.stackexchange.com/questions/261692/how-can-i-get-better-than-on2-space-complexity-for-a-type-sequence-search
		*/
	private:

		/* No. of elements in tuple used as stopping criterion */
		static constexpr size_t N = std::tuple_size_v<_Tuple>;

		template< size_t, bool, size_t... > struct recurse;

		/* SPECIALIZATION 1: Predicate is true, append to indices */
		/* Recursive case, non-empty list */
		template< size_t I, size_t... Is >
		struct recurse< I, true, Is... > {
			static constexpr bool pass = _Pred< std::tuple_element_t<I + 1, _Tuple > >::value;
			using type = typename recurse< I + 1, pass, Is..., I >::type;
		};
		/* Final case */
		template< size_t... Is >
		struct recurse< N - 1, true, Is... > {
			using type = std::index_sequence< Is..., N - 1 >;
		};

		/* SPECIALIZATION 2: Predicate is false, continue */
		/* Recursive case, non-empty list */
		template< size_t I, size_t... Is >
		struct recurse< I, false, Is... > {
			static constexpr bool pass = _Pred< std::tuple_element_t<I + 1, _Tuple > >::value;
			using type = typename recurse< I + 1, pass, Is... >::type;
		};
		/* Final case */
		template< size_t... Is >
		struct recurse< N - 1, false, Is... > {
			using type = std::index_sequence< Is... >;
		};

		static constexpr bool pass0 = _Pred< std::tuple_element_t<0, _Tuple > >::value;

	public:
		using type = typename recurse< 0, pass0 >::type;
	};

	/* SPECIALIZATION: EMPTY TUPLE */
	template <template<typename> typename _Pred>
	class indices_if< _Pred, std::tuple<> > {
	public:
		using type = std::index_sequence<>;
	};

/********************************************************************/
/* REPACK_IF: Build a variadic class from types satisfying a predicate */

	private:

	public:
		template< template<typename> typename _Predicate, typename...>
		struct repack_if;

		/* Non-empty tuple case */
		template< template<typename> typename _Pred, template<typename...> typename _VTC, typename T, typename... Ts>
		class repack_if< _Pred, _VTC<T, Ts...> > {
			using as_tuple_t = typename repack_as<std::tuple, _VTC<T, Ts...> >::type;
			using indices = typename indices_if<_Pred, as_tuple_t>::type;
			using sub_t = typename subtuple<as_tuple_t, indices>::type;
		public:
			using type = typename repack_as< _VTC, sub_t >::type;
		};

		/* Empty tuple specialization */
		template< template<typename> typename _Pred, template<typename...> typename _VTC>
		class repack_if< _Pred, _VTC<> > {
		public:
			using type = _VTC<>;
		};
	};
}