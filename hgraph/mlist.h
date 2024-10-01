#pragma once

#include <array>
#include <type_traits>
#include <memory>

/*****************************************************************************************/
/* INTRUSIVE MANAGED MULTILIST */

/*
MOTIVATION.
	Bidirectional iteration over arbitrary objects with always-valid iterators.

FEATURES.
	- multidimensionality: an element can be referred by any number of multilists.
	- iterator safety: default iterator is never invalidated by list operations 
	including insertion or removal.
	- STL-compliant range-for loop uses safe iterator by default.
	- behaves as a set, can only store unique elements
	- intrusiveness: element removal does not require lookup

BASIC USAGE.
	- stored elements must inherit mlist_hook publicly and declare their type using CRTP:
		class MyElement : public mlist_hook<MyElement> {...}
	- list declaration: list< MyElement >
	- standard list operations include push_back, push_front, erase etc
	- basic unsafe iterator provides ubegin(), uend(), urbegin(), urend()
	- basic unsafe iterator is *not* compatible with range_based for-loop, safe by default
	- iterator can be checked against nullptr

MULTIDIMENSIONAL HOOKS.
	- the same element can be referenced by N multilists
	- provide N in inherited hook declaration, eg:
		class MyElement : public mlist_hook<MyElement, 2> {...}
	- provide dimension I in multilist declaration:
		list<MyElement, 0> iterates the first dimension
		list<MyElement, 1> iterates the second dimension
	- I < N compatibility is checked at compile time

SAFE (MANAGED) ITERATION.
	- safe iteration relies on std::shared_ptr
	- stored elements inherit std::enable_shared_from_this through mlist_hook inheritance
	- stored elements must already be managed by a std::shared_ptr before access to
		a safe iterator
	- safe iterator provides default begin(), end(), rbegin(), rend()
	- range for-loop uses safe iterator by default

KNOWN QUIRKS.
	- using range-for loop without reference breaks shared_ptr communication,
	leading to double deletion and SEH exception. Always use
		for(auto& it : mlist),
	never
		for(auto it : mlist)

UNIMPLEMENTED. Features not used in hgraph that might be useful in general:
	- postfix ++/--, eg for element removal during unsafe iteration
	- insert_after and insert_before
*/

namespace hgraph {

	/* Tags for type safety */
	struct mlist_hook_tag {};
	struct mlist_tag {};

/************************************************************************************************/
/* INTRUSIVE HOOK */

		/* Basic hook for intrusive list. <_N> specifies multiple hooks (for use in edges and cliques) */
	template<typename T, size_t _N = 1>
	class mlist_hook :
		mlist_hook_tag,
		public std::enable_shared_from_this< T >
	{
		template<typename, size_t> friend class mlist;
	public:
		/* Number of dimensions */
		static constexpr size_t hook_dimensionality = _N;
		mlist_hook() : next{ nullptr }, prev{ nullptr } {};
	private:
		/* Array of pointers to next elements */
		std::array< T*, hook_dimensionality > next;
		/* Array of pointers to previous elements */
		std::array< T*, hook_dimensionality > prev;
	};

/************************************************************************************************/
/* INTRUSIVE MULTILIST */

		/* Intrusive (multi-)list. For multiple-hook elements, optional argument <_I> specifies
		the index of the hook associated with the list. */
	template<typename T, size_t _I = 0>
	class mlist : mlist_tag {
	public:
		static_assert(std::is_base_of_v< mlist_hook_tag, T >);
		static constexpr size_t I = _I;
		using type = T;
		static_assert(I < T::hook_dimensionality);

	protected:
		/* Pointer to first element */
		T* front_;
		/* Pointer to back_ element */
		T* back_;
		/* List size */
		size_t size_;

	public:
		mlist() : front_{ nullptr }, back_{ nullptr }, size_{ 0 } {};

		/* Getters */
		inline T& front() const { return *front_; }
		inline T& back() const { return *back_; }

		inline T* front_ptr() const { return front_; }
		inline T* back_ptr() const { return back_; }

		inline size_t size() const { return size_; }

/************************************************************************************************/
/* ITERATION */

		/* Iterator interface with begin, end and prefix ++/-- */
		// https://stackoverflow.com/questions/7758580/writing-your-own-stl-container/7759622#7759622
		// https://stackoverflow.com/questions/8164567/how-to-make-my-custom-type-to-work-with-range-based-for-loops
		// https://medium.com/@vgasparyan1995/how-to-write-an-stl-compatible-container-fc5b994462c6

	private:
		/* Base class, safe iterator based on shared_ptr*/
		class iterator_safe_base {
		public:
			using type = T;
			static constexpr size_t I = _I;
			bool operator!=(const iterator_safe_base& other) const {
				return ptr != other.ptr;
			}
			T& operator*() const { return *ptr; }
			T* operator->() const { return &*ptr; }
			inline bool is_null() { return ptr == nullptr; }

			/* Access some shared_ptr members */
			inline auto use_count() { return ptr.use_count(); }
			inline auto reset() { return ptr.reset(); }

		protected:
			iterator_safe_base(T* ptr) : ptr{ ptr ? ptr->shared_from_this() : nullptr } {}
			std::shared_ptr< T > ptr;
		};

	public:
		/* Forward safe iterator */
		class iterator_safe : public iterator_safe_base {
		public:
			iterator_safe(T* ptr) : iterator_safe_base(ptr) {}
			iterator_safe& operator++() {
				iterator_safe_base::ptr = iterator_safe_base::ptr->next[I] ?
					iterator_safe_base::ptr->next[I]->shared_from_this() : nullptr;
				return *this;
			}
			iterator_safe& operator--() {
				iterator_safe_base::ptr = iterator_safe_base::ptr->prev[I] ? 
					iterator_safe_base::ptr->prev[I]->shared_from_this() : nullptr;
				return *this;
			}
		};

		/* Reverse safe iterator */
		class reverse_iterator : public iterator_safe_base {
		public:
			reverse_iterator(T* ptr) : iterator_safe_base(ptr) {}
			reverse_iterator& operator--() {
				iterator_safe_base::ptr = iterator_safe_base::ptr->next[I] ?
					iterator_safe_base::ptr->next[I]->shared_from_this() : nullptr;
				return *this;
			}
			reverse_iterator& operator++() {
				iterator_safe_base::ptr = iterator_safe_base::ptr->prev[I] ? 
					iterator_safe_base::ptr->prev[I]->shared_from_this() : nullptr;
				return *this;
			}
		};

	public:
		/* Protected-range forward iterator */
		class iterator_ranged {
		public:
			using type = T;
			static constexpr size_t I = _I;
			bool operator!=(const iterator_safe_base& other) const {
				return ptr != other.ptr;
			}

			iterator_ranged(T* ptr, T* stop) :
				ptr{ ptr ? ptr->shared_from_this() : nullptr },
				prev{ ptr ? ptr->prev[I] : nullptr },
				stop{ stop ? stop->shared_from_this() : nullptr } {}

			/* Test whether current state matches stop criterion */
			explicit operator bool() const {				
				return (ptr != nullptr) && (stop.get() != prev);
			}

			iterator_ranged& operator++() {
				ptr = ptr->next[I] ? ptr->next[I]->shared_from_this() : nullptr;
				prev = ptr.get();
				return *this;
			}

			T& operator*() const { return *ptr; }
			T* operator->() const { return &*ptr; }
			inline bool is_null() { return ptr == nullptr; }
			/* Access some shared_ptr members */
			inline auto use_count() { return ptr.use_count(); }
			inline auto reset() { 
				stop.reset();
				return ptr.reset();
			}

		protected:
			/* Pointer to list element before current */
			T* prev;
			/* Managed ptr to current element */
			std::shared_ptr< T > ptr;
			/* Protected ptr to stopping element */
			const std::shared_ptr< T > stop;
		};

	private:
		/* Base class, unsafe iterator based on raw ptr*/
		class iterator_unsafe_base {
		public:
			using type = T;
			static constexpr size_t I = _I;
			iterator_unsafe_base(T* ptr) : ptr{ ptr } {}
			bool operator!=(const iterator_unsafe_base& other) const {
				return ptr != other.ptr;
			}
			T& operator*() const { return *ptr; }
			T* operator->() const { return ptr; }
			inline bool is_null() { return ptr == nullptr; }
		protected:
			T* ptr;
		};

	public:
		/* Forward unsafe iterator */
		class iterator_unsafe : public iterator_unsafe_base {
		public:
			iterator_unsafe(T* ptr) : iterator_unsafe_base(ptr) {}
			iterator_unsafe& operator++() {
				iterator_unsafe_base::ptr = iterator_unsafe_base::ptr->next[I];
				return *this;
			}
			iterator_unsafe& operator--() {
				iterator_unsafe_base::ptr = iterator_unsafe_base::ptr->prev[I];
				return *this;
			}
		};

		/* Reverse unsafe iterator */
		class reverse_iterator_unsafe : public iterator_unsafe_base {
		public:
			reverse_iterator_unsafe(T* ptr) : iterator_unsafe_base(ptr) {}
			reverse_iterator_unsafe& operator--() {
				iterator_unsafe_base::ptr = iterator_unsafe_base::ptr->next[I];
				return *this;
			}
			reverse_iterator_unsafe& operator++() {
				iterator_unsafe_base::ptr = iterator_unsafe_base::ptr->prev[I];
				return *this;
			}
		};

	public:

		/* Safe (default) iterator */
		iterator_safe begin() const { return iterator_safe(front_); }
		iterator_safe end()   const { return iterator_safe(nullptr); }
		iterator_safe last() const { return iterator_safe(back_); }

		reverse_iterator rbegin() const { return reverse_iterator(back_); }
		reverse_iterator rend()   const { return reverse_iterator(nullptr); }

		/* Range-protected iterator */
		iterator_ranged gbegin() const { return iterator_ranged(front_, back_); }

		/* Unsafe iterator */
		iterator_unsafe ubegin() const { return iterator_unsafe(front_); }
		iterator_unsafe uend()   const { return iterator_unsafe(nullptr); }

		reverse_iterator_unsafe urbegin() const { return reverse_iterator_unsafe(back_); }
		reverse_iterator_unsafe urend()   const { return reverse_iterator_unsafe(nullptr); }

/************************************************************************************************/
/* ELEMENT INSERTION */

		/* Append element to list */
		inline void push_back(T* x) {
			static_assert(std::is_base_of_v<mlist_hook_tag, T>);
			assert(x != nullptr);

			/* Attempt to reinsert an element ? */
			assert( (x->prev[I] == nullptr) && (x->next[I] == nullptr) &&
				"Attempt to insert a duplicate element."
			);

			/* First element in list ? */
			if (front_ == nullptr) {
				assert(back_ == nullptr);
				front_ = x;
				back_ = x;
			}
			else {
				assert(back_ != nullptr);
				back_->next[I] = x;
				x->prev[I] = back_;
				back_ = x;
			}

			size_++;
		}

		/* Prepend element to list */
		inline void push_front(T* x) {
			static_assert(std::is_base_of_v<mlist_hook_tag, T>);
			assert(x != nullptr);

			/* Attempt to reinsert an element ? */
			assert((x->prev[I] == nullptr) && (x->next[I] == nullptr) &&
				"Attempt to insert a duplicate element."
			);

			/* First element in list ? */
			if (front_ == nullptr) {
				assert(back_ == nullptr);
				front_ = x;
				back_ = x;
			}
			else {
				assert(back_ != nullptr);
				front_->prev[I] = x;
				x->next[I] = front_;
				front_ = x;
			}

			size_++;
		}

/************************************************************************************************/
/* ELEMENT REMOVAL */

		/* Remove element from list */
		inline void erase(T* x) {
			static_assert(std::is_base_of_v<mlist_hook_tag, T>);
			assert(x != nullptr);

			if (front_ == x) front_ = x->next[I];
			if (back_ == x) back_ = x->prev[I];

			if (x->prev[I] != nullptr) x->prev[I]->next[I] = x->next[I];
			if (x->next[I] != nullptr) x->next[I]->prev[I] = x->prev[I];

			/* Cleanup hooks to enable reinsertion */
			x->prev[I] = nullptr;
			x->next[I] = nullptr;

			size_--;
		}

		/* Remove element from list */
		inline void erase(T& x) {
			static_assert(std::is_base_of_v<mlist_hook_tag, T>);
			erase(&x);
		}

		/* Remove element from list */
		inline void erase(iterator_safe& x) {
			static_assert(std::is_base_of_v<mlist_hook_tag, T>);
			erase(&*x);
		}

		/* Remove element from list */
		inline void erase(iterator_unsafe& x) {
			static_assert(std::is_base_of_v<mlist_hook_tag, T>);
			erase(&*x);
		}

/************************************************************************************************/
/* LIST REORDERING */

		/* Return a new vector of pointers to list elements */
		inline std::vector<T*> as_vector() {

			std::vector<T*> v;
			for(auto it = ubegin(); it != uend(); ++it) {
				v.push_back( &*it );
			}

			return v;
		}

		/* Reorder list along a vector of pointers to its elements.
		Typical usage is generate vector*/
		inline void reorder(const std::vector<T*> v) {

			assert(v.size() == size());

			if (v.size() == 0) return;
			
			for (auto it : v) {
				erase(it);
			}

			for (auto it : v) {
				push_back(it);
			}
		}

	}; // class mlist
} // namespace hgraph