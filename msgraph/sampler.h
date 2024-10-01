#pragma once

#include <vector>
#include <unordered_map>

class Sampler {
public:
	void seed(uint32_t s1, uint32_t s2);

	Sampler() {
		seed(123, 456);
	}

	Sampler(uint32_t s1, uint32_t s2) {
		seed(s1, s2);
	};

	/* HACK hash functor for vectors, can do better ;-)*/
	struct vector_hash {
		size_t operator()(const std::vector< size_t >& v) const {
			size_t h = 0;
			for (size_t i = 0; i < v.size(); ++i) {
				h ^= std::hash<size_t>()(v[i] << 1);
			}
			return h;
		}
	};

	/* Binomial random deviate */
	size_t rbinom(const size_t n, const double p);

	/* Multinomial random deviate */
	std::vector< size_t > rmultinom(const size_t n, std::vector< double >& p);

	/* Multivariate hypergeometric random deviate */
	std::vector< size_t > rmvhyper(const size_t n, std::vector< size_t >& K);

	/* Multivariate hypergoemtric probability mass */
	static double dmvhyper(std::vector< size_t > k, const size_t n, std::vector< size_t > K, const size_t N);

	/* Multidimensional hypergeometric random deviate (full S x T contingency table with
	T columns, constant column sum n and variable row sums K) */
	std::vector< std::vector< size_t > > rctable(const size_t n, const size_t T, std::vector< size_t >& K);
	
	/* Approximate multidimensional hypergeometric random deviate */
	std::unordered_map< std::vector< size_t >, size_t, vector_hash >
		rctable_approx(const size_t n, const size_t T, std::vector< size_t >& K, const double precision = 1.);
};