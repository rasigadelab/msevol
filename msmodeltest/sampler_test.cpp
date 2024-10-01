#include "pch.h"

#include <vector>
//#include <algorithm>
#include <numeric>
#include <iostream>
#include <unordered_map>

using std::cout;
using std::endl;

//https://github.com/wch/r-source/tree/trunk/src/nmath

#define MATHLIB_STANDALONE 1
#include "../vendor/rmath/Rmath.h"

#define normInv(x) (qnorm((x), 0.0, 1.0, 1, 0))
#define normRand() (rnorm(0.0, 1.0))

/* Unified sampler class */
/* 
Should implement binomial sampling, hypergeometric sampling,
both exact and with tunable approximations. Note that rmath routines 
rbinom and rhyper take and return double.
*/

#include "../msgraph/sampler.h"

TEST(Sampler, rbinom) {
	Sampler sp;

	const double p = 0.2;
	const size_t n = 5;

	for (size_t i = 0; i < 10; ++i) {
		size_t k = sp.rbinom(n, p);
		//printf("Draw result = %zi\n", k);
	}
}

TEST(Sampler, rmultinom) {
	Sampler sp;

	std::vector<double> p = { 0.2, 0.3, 0.5 };

	const size_t n = 5;

	for (size_t i = 0; i < 10; ++i) {

		std::vector< size_t > draw = sp.rmultinom(n, p);
		ASSERT_EQ(std::accumulate(draw.cbegin(), draw.cend(), size_t(0)), n);
		//cout << "Draw:" << draw[0] << draw[1] << draw[2] << endl;
		//printf("Draw result = %zi, %zi, %zi\n", draw[0], draw[1], draw[2]);

	}
}

TEST(Sampler, rmvhyper) {

	Sampler sp;

	std::vector< size_t > v = { 10, 3, 4 }; // Population

	size_t n = 5; // No. of draws

	for (size_t i = 0; i < 10; ++i) {

		std::vector< size_t > draw = sp.rmvhyper(n, v);
		ASSERT_EQ(std::accumulate(draw.cbegin(), draw.cend(), size_t(0)), n);
		//cout << "Draw:" << draw[0] << draw[1] << draw[2] << endl;
		//printf("Draw result = %zi, %zi, %zi\n", draw[0], draw[1], draw[2]);

	}
}

TEST(Sampler, rctable) {

	/* Dispatch hypergeometric distribution */
	/* This is the distribution obtained by dispatching a population of N balls with S categories
	into T bins of equal size n, with n*T = S, until all balls have been drawn */

	Sampler sp;

	std::vector< size_t > K = { 10, 30, 40}; // Population
	size_t S = K.size();
	size_t N = std::accumulate(K.cbegin(), K.cend(), size_t(0));

	size_t T = 5; // No. of bins
	size_t n = 16;

	std::vector< std::vector< size_t > > v = sp.rctable(n, T, K);

	/* Check margins */
	std::vector< size_t > rowSums(S);
	std::vector< size_t > colSums(T);
	for (size_t j = 0; j < T; ++j) {
		for (size_t i = 0; i < S; ++i) {
			rowSums[i] += v[j][i];
			colSums[j] += v[j][i];
		}
	}

	for (size_t j = 0; j < T; ++j) {
		ASSERT_EQ(colSums[j], n);
	}

	for (size_t i = 0; i < S; ++i) {
		ASSERT_EQ(rowSums[i], K[i]);
	}
}

TEST(Sampler, dmvhyper) {

	std::vector< size_t > K = { 100, 300, 400 }; // Population
	std::vector< size_t > k = { 5, 11, 13 }; // Draw
	size_t N = std::accumulate(K.cbegin(), K.cend(), size_t(0)); // Population size
	size_t n = std::accumulate(k.cbegin(), k.cend(), size_t(0)); // Draw size

	double p = Sampler::dmvhyper(k, n, K, N);

	//printf("dmvhyper prob = %.3e\n", p);

	/* Check against R result */
	ASSERT_DOUBLE_EQ(std::round(p * 1e5), std::round(0.023455758631108267 * 1e5));
}

TEST(Sampler, rctable_approx) {

	Sampler sp;
	sp.seed(111, 123);

	std::vector< size_t > K = { 796, 2, 2 }; // Population
	size_t S = K.size();
	size_t T = 50; // No. of bins
	size_t n = 16;

	auto res = sp.rctable_approx(n, T, K, 2.);

	/* Print results */
	for (auto& it : res) {
		//cout << "Draw = " << it.first[0] << " " << it.first[1] << " " << it.first[2] << " - h = " << it.second << endl;
	}

	/* Sum of hashmap's h should equal T */
	size_t h_check = 0;
	for (auto& it : res) {
		h_check += it.second;
	}
	EXPECT_EQ(h_check, T);

	/* All categories must be dispatched exactly */
	for (size_t i = 0; i < S; ++i) {
		size_t k_check = 0;
		for (auto& it : res) {
			k_check += it.first[i] * it.second;
		}
		EXPECT_EQ(k_check, K[i]);
	}

	/* Remark that algorithm stops very fast if there are few events to be consumed.
	No need to add any check on that. */
}