#include "pch.h"
#include "../hgraph/mlist.h"
#include <memory>
#include <random>
#include <algorithm>

using namespace hgraph;

class Hooked : public mlist_hook<Hooked> {
public:
	const size_t id;
	Hooked(size_t id_ = 0) : id{ id_ } {}
};
class MultiHooked : public mlist_hook<MultiHooked, 2> {
public:
	const size_t id;
	MultiHooked(size_t id_ = 0) : id{ id_ } {}
};
class Root {
public:
	mlist<Hooked> list;
	mlist<MultiHooked, 0> mlist0;
	mlist<MultiHooked, 1> mlist1;
};

TEST(mlist, push_back) {

	Root r;

	ASSERT_EQ(r.list.front_ptr(), nullptr);
	ASSERT_EQ(r.list.back_ptr(), nullptr);
	ASSERT_EQ(r.list.size(), 0);

	ASSERT_EQ(r.mlist0.front_ptr(), nullptr);
	ASSERT_EQ(r.mlist0.back_ptr(), nullptr);
	ASSERT_EQ(r.mlist0.size(), 0);

	ASSERT_EQ(r.mlist1.front_ptr(), nullptr);
	ASSERT_EQ(r.mlist1.back_ptr(), nullptr);
	ASSERT_EQ(r.mlist1.size(), 0);

	const size_t n = 5;
	for (size_t i = 0; i < n; i++) {
		r.list.push_back(new Hooked(i));
	}
	ASSERT_EQ(r.list.size(), n);

	/* Element order should be preserved */
	size_t expected_id = 0;
	for(auto it = r.list.ubegin(); it != r.list.uend(); ++it) {
		ASSERT_EQ(it->id, expected_id++);
	}
	ASSERT_EQ(expected_id, n);

	/* Reverse unsafe iterator */
	expected_id = n;
	for(auto it = r.list.urbegin(); it != r.list.urend(); ++it) {
		ASSERT_EQ(it->id, --expected_id);
	}
	ASSERT_EQ(expected_id, 0);
}

TEST(mlist, erase) {

	Root r;

	const size_t n = 5;
	for (size_t i = 0; i < n; i++) {
		r.list.push_back(new Hooked(i));
	}

	/* Remove element #2 */
	auto h = r.list.ubegin();
	++h;
	++h;
	ASSERT_EQ(h->id, 2);

	/* Erase by raw pointer */
	r.list.erase(h);

	for (auto it = r.list.ubegin(); it != r.list.uend(); ++it) {
		ASSERT_NE(it->id, 2);
	}
}

TEST(mlist, iteratorUnsafe) {
	Root r;

	/* Single-hook list */
	const size_t n = 5;
	for (size_t i = 0; i < n; i++) {
		r.list.push_back(new Hooked(i));
	}

	/* Unsafe iteration */
	size_t i_check = 0;
	for (auto it = r.list.ubegin(); it != r.list.uend(); ++it) {
		ASSERT_EQ(it->id, i_check++);
	}
	ASSERT_EQ(i_check, r.list.size());

	/* Multi-hook list, fill each hook by opposite order */
	for (size_t i = 0; i < n; i++) {
		MultiHooked* m = new MultiHooked(i);
		r.mlist0.push_back(m);
		r.mlist1.push_front(m);
	}

	/* List 0 should be filled in forward order */
	i_check = 0;
	for (auto it = r.mlist0.ubegin(); it != r.mlist0.uend(); ++it) {
		ASSERT_EQ(it->id, i_check++);
	}
	ASSERT_EQ(i_check, r.mlist0.size());

	/* List 1 should be filled in reverse order */
	i_check = r.mlist1.size();
	for (auto it = r.mlist1.ubegin(); it != r.mlist1.uend(); ++it) {
		ASSERT_EQ(it->id, --i_check);
	}
	ASSERT_EQ(i_check, 0);
}

TEST(mlist, iteratorSafe) {
	
	/* Single-hook list */
	Root r;
	const size_t n = 5;

	/* Managed shared_ptr vector, required for safe list iteration */
	std::vector< std::shared_ptr< Hooked > > v;
	v.reserve(n);

	for (size_t i = 0; i < n; i++) {
		v.push_back(std::make_shared< Hooked >(i));
		r.list.push_back( v[i].get() );
	}

	/* Forward iteration, range-based for loop */
	size_t i_check = 0;
	for (auto& it : r.list) {
		ASSERT_EQ(it.id, i_check++);
	}

	/* Forward iteration, iterator-based */
	i_check = 0;
	for (auto it = r.list.begin(); it != r.list.end(); ++it) {
		ASSERT_EQ(it->id, i_check++);
	}
	ASSERT_EQ(i_check, n);

	/* Reverse iteration, iterator-based */
	i_check = n;
	for (auto it = r.list.rbegin(); it != r.list.rend(); ++it) {
		ASSERT_EQ(it->id, --i_check);
	}
	ASSERT_EQ(i_check, 0);
}

// WIP
TEST(mlist, iterateProtectedRange) {
	/* Protected range iteration guarantees that 
	visited elements are only those that were present at the beginning
	of the iteration */
	
	/* Single-hook list */
	Root r;
	const size_t n = 5;

	/* Managed shared_ptr vector, required for safe list iteration */
	std::vector< std::shared_ptr< Hooked > > v;
	v.reserve(n);

	for (size_t i = 0; i < n; i++) {
		v.push_back(std::make_shared< Hooked >(i));
		r.list.push_back(v[i].get());
	}

	/* Protected-range iteration controls against a protected stopper that points to
	the last element in list. */
	size_t i_check = 0;
	for (auto [it, stop] = std::tuple{ r.list.begin(), r.list.last() }; it != stop; ++it) {
		ASSERT_EQ(it->id, i_check++);
	}
	ASSERT_EQ(i_check, n - 1);

	/* Protected-range iterator tracks both current and previous element + 
	locks the test element to prevent its deletion + checks equality against prev
	pointer */

	i_check = 0;
	for (auto it = r.list.gbegin(); it; ++it) {
		ASSERT_EQ(it->id, i_check++);
	}
	ASSERT_EQ(i_check, n - 1);

	/* Range-protected iteration should stop even if list grows */
	i_check = 0;
	for (auto it = r.list.gbegin(); it; ++it) {
		/* Add a new element */
		v.push_back(std::make_shared< Hooked >(n + i_check));
		r.list.push_back(v[n + i_check].get());
		ASSERT_EQ(it->id, i_check++);
	}
	ASSERT_EQ(i_check, n - 1);
	ASSERT_EQ(r.list.size(), 2*n - 1);



}

TEST(mlist, iteratorKillstop) {
/* Protected range iteration guarantees that
visited elements are only those that were present at the beginning
of the iteration */

/* Single-hook list */
	Root r;
	const size_t n = 5;

	/* Managed shared_ptr vector, required for safe list iteration */
	std::vector< std::shared_ptr< Hooked > > v;
	v.reserve(n);

	for (size_t i = 0; i < n; i++) {
		v.push_back(std::make_shared< Hooked >(i));
		r.list.push_back(v[i].get());
	}

	/* WIP: Range-protected iteration should stop even if stop element is removed for owner */
	auto killstop = [&v]() {
		v.erase(v.begin() + v.size() - 1);
	};
	size_t i_check = 0;
	for (auto it = r.list.gbegin(); it; ++it) {
		if (i_check == 0) killstop();
		/* Add a new element */
		v.push_back(std::make_shared< Hooked >(n + i_check - 1));
		r.list.push_back(v[n + i_check - 1].get());
		ASSERT_EQ(it->id, i_check++);
	}
	ASSERT_EQ(i_check, n - 1);
	ASSERT_EQ(r.list.size(), 2 * n - 1);
}

TEST(mlist, shuffle) {

	/* Randomly shuffle all elements in the list
	* Naive implementation stores order as vector of pointers, shuffles that vector then rebuilds list
	*/

	Root r;
	const size_t n = 5;

	/* Managed shared_ptr vector, required for safe list iteration */
	std::vector< std::shared_ptr< Hooked > > v;
	v.reserve(n);

	for (size_t i = 0; i < n; i++) {
		v.push_back(std::make_shared< Hooked >(i));
		r.list.push_back(v[i].get());
	}

	/* Build hook vector (raw pointers) and shuffle it, then reorder list */

	auto rv = r.list.as_vector();

	/* Shuffle vector of pointers */
	std::random_device rd;
	std::mt19937 g(rd());
	std::shuffle(rv.begin(), rv.end(), g);

	/* Reorder list along vector of pointers */
	r.list.reorder(rv);

	ASSERT_EQ(r.list.size(), n);

	/* Vector should reflect list order*/
	auto rvit = rv.begin();	
	for (auto it = r.list.begin(); it != r.list.end(); ++it) {
 		//std::cout << it->id  << " " << (*rvit++)->id << "  ";
		ASSERT_EQ(it->id, (*rvit++)->id);
	}
	//std::cout << std::endl;
}