#pragma once

#include <tuple>
#include <type_traits>

/* Metaprogramming class */

namespace hgraph {

	/* Minimal placeholder */
	template<typename...>
	class tl {};

	template<template<typename...> typename, typename>
	class typelist;

	template<template<typename...> typename _VTC, typename... _Types>
	class typelist< _VTC, _VTC< _Types... > > {
	public:
		/* Direct access to type */
		using type = _VTC< _Types... >;

		/* Length of parameter pack */
		static constexpr size_t size = sizeof...(_Types);

		/* Apply a metafunction on each type */
		template< template<typename> typename _MetaFunc >
		struct transform {
			using type = _VTC< typename _MetaFunc< _Types >::type... >;
			using chain = typelist< _VTC, type >;
		};

		/* Change holder type */
		template< template<typename...> typename _NewVTC >
		struct repack_as {
			using type = _NewVTC< _Types... >;
			using chain = typelist< _NewVTC, type >;
		};

		/* Conditional, range-protected access to types */
	protected:
		template< size_t, bool > struct at_conditional;
		template< size_t I >
		struct at_conditional< I, true > {
			using type = std::tuple_element_t< I, std::tuple< _Types... > >;
		};
		template< size_t I >
		struct at_conditional< I, false > {
			using type = void;
		};

	public:
		/* Get Nth type, return void is out of bounds */
		template<size_t I>
		struct at : at_conditional < I, I < size > {};

		/* Append to list */
		template< typename T >
		struct push_back {
			using type = _VTC< _Types..., T>;
			using chain = typelist< _VTC, type >;
		};

		/* Prepend to list */
		template< typename T >
		struct push_front {
			using type = _VTC< T, _Types...>;
			using chain = typelist< _VTC, type >;
		};

	private:
		/* Packs bool to help expansion */
		template<bool _B>
		struct expand_bool {
			static constexpr bool value = _B;
		};


	public:
		/* Extract types based on integer sequence */
		template< typename > class extract;
		template< size_t... Is >
		class extract< std::index_sequence< Is... > > {
			static_assert(std::conjunction_v < expand_bool< Is < size >... >,
				"Typelist index out of bounds in sequence.");
		private:
			using as_tuple = std::tuple< _Types... >;
		public:
			using type = _VTC< typename std::tuple_element<Is, as_tuple>::type... >;
			using chain = typelist< _VTC, type >;
		};

		/* Find indices of types satifying a predicate */
		template< template<typename> typename _Pred >
		class indices_of {
		private:
			using as_tuple = std::tuple< _Types... >;

			template< size_t, bool, size_t... > struct recurse;

			/* SPECIALIZATION 1: Predicate is true, append to indices */
			/* Recursive case, non-empty list */
			template< size_t I, size_t... Is >
			struct recurse< I, true, Is... > {
				static constexpr bool pass = _Pred< std::tuple_element_t<I + 1, as_tuple > >::value;
				using type = typename recurse< I + 1, pass, Is..., I >::type;
			};
			/* Final case */
			template< size_t... Is >
			struct recurse< size - 1, true, Is... > {
				using type = std::index_sequence< Is..., size - 1 >;
			};

			/* SPECIALIZATION 2: Predicate is false, continue */
			/* Recursive case, non-empty list */
			template< size_t I, size_t... Is >
			struct recurse< I, false, Is... > {
				static constexpr bool pass = _Pred< std::tuple_element_t<I + 1, as_tuple > >::value;
				using type = typename recurse< I + 1, pass, Is... >::type;
			};
			/* Final case */
			template< size_t... Is >
			struct recurse< size - 1, false, Is... > {
				using type = std::index_sequence< Is... >;
			};

		public:
			static constexpr bool pass = _Pred< std::tuple_element_t<0, as_tuple > >::value;
			// Fixme, treat empty case
			//using type = typename recurse< 0, pass >::type;

			using type = typename std::enable_if<
				(size > 0),
				typename recurse< 0, pass >::type
			>::type;
		};

		/* Filter types that satisfy a predicate */
		template< template<typename> typename _Pred >
		class filter : public extract< typename indices_of< _Pred >::type > {};

		/* Check whether types contain a given type */
		template< typename T >
		struct has_type : std::disjunction< std::is_same< T, _Types >... > {};

	};

	/* Wrapper for std::tuple */
	template< typename > class typetuple;
	template< typename... Ts >
	class typetuple< std::tuple< Ts... > > : public typelist< std::tuple, std::tuple< Ts... > > {};

}