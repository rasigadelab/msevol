#pragma once

#include <array>
#include <string>
#include <tuple>
using namespace std;

#include "TableFile.h"

class Writer : public TableFile {
private:
	/* Helpers */
	inline void newSeparator() {
		fprintf(file, separator);
	}
	inline void newLine() {
		fprintf(file, "\n");
	}

	/* Serialize */
	template<typename T1, typename T2, typename... T>
	inline void serialize(T1& currentValue, T2& nextValue, T&... rest) {
		write(currentValue);
		newSeparator();
		return serialize(nextValue, rest...);
	};
	template<typename T1>
	inline void serialize(T1& currentValue) {
		write(currentValue);
	};

	/* Individual writers */
	inline void write(const size_t& x) {
		fprintf(file, "%Iu", x);
	}
	inline void write(const double& x) {
		fprintf(file, "%.16G", x);
	}
	inline void write(const char* x) {
		fprintf(file, "%s", x);
	}
	inline void write(const string& x) {
		fprintf(file, "%s", x.c_str());
	}
	/* Array writer */
	template<size_t I, typename T, size_t N>
	inline void write_array(const array<T, N>& x) {
		write(x[I]);
		if constexpr ((I + 1) < N) {
			newSeparator();
			write_array<I + 1>(x);
		}
	}
	template<typename T, size_t N>
	inline void write(const array<T, N>& x) {
		if constexpr (N > 0) {
			write_array<0>(x);
		}
	}

protected:
	/* Specify file name */
	Writer(const char* fname_in, const bool append = false)
		: TableFile(fname_in, append ? "ab" : "wb") {}

	/* Automatic file name */
	Writer(const char* path, const char* fname_in, const bool append = false)
		: TableFile(path, fname_in, append ? "ab" : "wb") {}

	~Writer() {}

	/* Set first line */
	template<typename... T>
	inline void firstline(T... header) {
		serialize(header...);
	}

	/* Tuple version */
	template<typename... T>
	inline void firstline_tuple(std::tuple< T... > header_tuple) {
		std::apply([this](auto&... h) { this->firstline(h...); }, header_tuple);
	}

	/* Append a line */
	template<typename... T>
	inline void push(T&... header) {
		newLine();
		serialize(header...);
	}

	/* Tuple version */
	template<typename... T>
	inline void push_tuple(std::tuple< T... > line_tuple) {
		std::apply([this](auto&... h) { this->push(h...); }, line_tuple);
	}
};
