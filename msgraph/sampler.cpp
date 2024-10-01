#include <numeric>
#include <cassert>
#include <algorithm>
#include <stdexcept>

#include "sampler.h"


//https://github.com/wch/r-source/tree/trunk/src/nmath

#define MATHLIB_STANDALONE 1
#include "../vendor/rmath/Rmath.h"

/* Disambiguate rbinom */
static inline double rmath_rbinom(double n, double p) {
	return rbinom(n, p);
}



/*
Numerically stable implementation of multivariate hypergeometric probability mass

for n draws of population K with S categories with count K_1, ..., K_S
and values k_1, ..., k_s

let N = sum(K)
let p = n / N

probability is

dbinom(k_1, K_1, p) * dbinom(k_2, K_2, p) * ... * dbinom(k_S, K_S, p) / dbinom(n, N, p)

which is probability better computed in log space given the amount of powers required.

signature of dbinom
double dbinom(double x, double n, double p, int give_log)

*/
double Sampler::dmvhyper(std::vector< size_t > k, const size_t n, std::vector< size_t > K, const size_t N) {

	size_t S = k.size();
	assert(S == K.size());
	double p = double(n) / double(N);

	double prob = 1. / dbinom(double(n), double(N), p, false);
	for (size_t i = 0; i < S; ++i) {
		prob *= dbinom(double(k[i]), double(K[i]), p, false);
	}
	/* Check for underflow */
	assert(p > 0.);
	return prob;
}

void Sampler::seed(uint32_t s1, uint32_t s2) {
	set_seed(s1, s2);
}

/* Binomial random deviate */
size_t Sampler::rbinom(const size_t n, const double p) {
	return size_t(rmath_rbinom(double(n), p));
}

/* Multinomial random deviate */
std::vector< size_t > Sampler::rmultinom(const size_t n, std::vector< double >& p) {

	const size_t s = p.size();
	std::vector< size_t > res(s);

	double ni = double(n);
	double pi = std::accumulate(p.cbegin(), p.cend(), 0.);

	for (size_t i = 0; i < s; ++i) {
		const double k = rmath_rbinom(ni, p[i] / pi);
		res[i] = size_t(std::round(k));
		ni -= k;
		pi -= p[i];
	}

	return res;
}

/* Multivariate hypergeometric random deviate */
std::vector< size_t > Sampler::rmvhyper(const size_t n, std::vector< size_t >& K) {

	const size_t s = K.size();	
	std::vector< size_t > res(s);
	size_t N = std::accumulate(K.cbegin(), K.cend(), size_t(0));
	assert(n <= N);

	if (s == 1) {
		res[0] = n;
		return res;
	}
	
	if (N == 0) {
		res = K;
		return res;
	}

	size_t ni = n;
	for (size_t i = 0; i < s; ++i) {
		const size_t k = size_t(std::round(rhyper(double(K[i]), double(N - K[i]), double(ni))));
		res[i] = k;
		N -= K[i]; // Remove from population
		ni -= k; // Remove from draws
	}

	return res;
}

/* Multidimensional hypergeometric random deviate (full S x T contingency table with
T columns, constant column sum n and variable row sums K) */
std::vector< std::vector< size_t > > Sampler::rctable(const size_t n, const size_t T, std::vector< size_t >& K) {

	const size_t S = K.size();

	/* Consistency check: population size must equal T*n, which results from
	msgraph design where no. of events to be dispatched is necessarily a multiple of
	container count (T) times container multiplicity (n) */
	const size_t N = std::accumulate(K.cbegin(), K.cend(), size_t(0));
	if (N != T * n) throw std::exception("Size mismatch in rctable.");

	std::vector< std::vector< size_t > > res(T);

	std::vector< size_t > Ki = K; /* Keep track of remaining balls */

	for (size_t j = 0; j < T; ++j) {
		res[j] = rmvhyper(n, Ki);
		/* Remove draws from population */
		auto& vj = res[j];
		for (size_t i = 0; i < S; ++i) {
			Ki[i] -= vj[i];
		}
	}

	return res;
}

/* Approximate multidimensional hypergeometric random deviate */
std::unordered_map< std::vector< size_t >, size_t, Sampler::vector_hash >
Sampler::rctable_approx(const size_t n, const size_t T, std::vector< size_t >& K, const double precision) {

	/* Precision factor: set > 1 for increased precision by reducing the probability of repeated draws;
	set < 1 for increased speed by increasing repeated draws */

	size_t S = K.size();
	size_t N = std::accumulate(K.cbegin(), K.cend(), size_t(0));
	if (N != T * n) throw std::exception("Size mismatch in rctable.");

	/* Result is a hashmap of vectors where value is the count. This is similar to barcode hashing */
	std::unordered_map< std::vector< size_t >, size_t, Sampler::vector_hash > res;

	std::vector< size_t > Ki = K; /* Keep track of remaining balls */
	size_t Ni = N; /* Keep track of remaining population */
	size_t Ti = T; /* Keep track of remaining columns (containers) */

	while (Ti > 0) {
		/* Draw a column */
		std::vector< size_t > k = rmvhyper(n, Ki);
		/* Probability of draw */
		double p = dmvhyper(k, n, Ki, Ni);

		/* Apply approximation factor */
		p = std::min(1., p / precision);

		/* Expected no. of identical draws */
		size_t h = rbinom(Ti, p);

		/* No. of draws cannot be zero if k was actually drawn */
		if (h == 0) h = 1;

		/* For each category, the h repeated draws cannot consume
		more than the remainder in Ki */
		for (size_t i = 0; i < S; ++i) {
			const size_t ki = k[i];
			if (ki == 0) continue;
			const size_t h_max = Ki[i] / ki;
			if (h > h_max) h = h_max;
		}

		/* Update counts in each category Ki, for total population Ni and no. of columns Ti */
		for (size_t i = 0; i < S; ++i) {
			Ki[i] -= k[i] * h;
		}
		Ni -= n * h;
		Ti -= h;

		/* Register draw in hashmap */
		auto it = res.find(k);
		if (it == res.end()) {
			/* New combination */
			res.insert({ std::move(k), h });
		}
		else {
			//cout << "Repeated combination !" << endl;
			/* Existing combination */
			it->second += h;
		}
	}
	assert(Ti == 0);
	assert(Ni == 0);

	return std::move(res);
}