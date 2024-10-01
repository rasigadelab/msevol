#pragma once

#include <array>
#include <string>
#include <tuple>
using namespace std;

#include "TableFile.h"
// https://docs.microsoft.com/en-us/cpp/c-runtime-library/format-specification-fields-scanf-and-wscanf-functions?view=vs-2019
// https://docs.microsoft.com/en-us/cpp/c-runtime-library/scanf-width-specification?view=vs-2019


class Reader : public TableFile {
protected:
	/**************************************************/
	template<typename T>
	struct Format {
		static char* read(T& x, char* ptr) {

			size_t offset = 0;

			if constexpr (std::is_same_v<T, double>) {
				x = stod(ptr, &offset);
			}
			else if constexpr (std::is_same_v<T, size_t>) {
				x = stoull(ptr, &offset);
			}
			else if constexpr (std::is_same_v<T, int>) {
				x = stoi(ptr, &offset);
			}
			else {
				static_assert(false, "Unimplemented format reader");
			}

			ptr += offset;

			/* Skip separator or blank */
			while (std::isblank(char(*ptr))) ptr++;

			return ptr;
		}
	};

	template< typename T, size_t N >
	struct Format< std::array<T, N> > {
		static char* read(std::array<T, N>& x, char* ptr) {
			for (size_t i = 0; i < N; i++) {

				ptr = Format<T>::read(x[i], ptr);

			}
			return ptr;
		}
	};

	/**************************************************/

	/* Line buffer */
	char* buffer;

	/* Buffer length */
	const size_t bufferLength;

	inline void init() {
		buffer = new char[bufferLength];
		fgets(buffer, (int)bufferLength, file);
	}

	/* Specify file name */
	Reader(const char* fname_in, size_t bufferLength_in = 65536)
		: TableFile(fname_in, "rb"), bufferLength{ bufferLength_in } {
		init();
	}
	/* Automatic file name */
	Reader(const char* path, const char* fname_in, size_t bufferLength_in = 65536)
		: TableFile(path, fname_in, "rb"), bufferLength{ bufferLength_in }  {
		init();
	}

	~Reader() {
		delete buffer;
	}

public:

	/* Tuple based alternative */
	template<size_t I, typename T>
	void read_tuple(T& t, char* ptr) {
		constexpr size_t N = std::tuple_size_v<T>;

		ptr = Format< std::tuple_element_t<I, T> >::read(std::get<I>(t), ptr);

		if constexpr ((I + 1) < N) return read_tuple< I + 1 >(t, ptr);
	}

	template<typename T>
	bool read_tuple(T& t) {
		char* ptr = fgets(buffer, (int)bufferLength, file);
		if (ptr == nullptr) return false;
		read_tuple<0>(t, buffer);
		return true;
	}

	/* Tuple based alternative with separate ID parsing */
	template<size_t I, size_t N, typename T>
	void read_tuple(T& t, char* ptr) {		
		ptr = Format< std::tuple_element_t<I, T> >::read(std::get<I>(t), ptr);
		if constexpr ((I + 1) < N) return read_tuple< I + 1, N >(t, ptr);
	}

	/* Vertex overload: read vertex id first */
	template<typename T>
	bool read_tuple(size_t& id, T& t) {
		/* Size of tuple, can be zero */
		constexpr size_t N = std::tuple_size_v<T>;
		char* ptr = fgets(buffer, (int)bufferLength, file);
		if (ptr == nullptr) return false;
		/* Read vertex ID */
		ptr = Format< size_t >::read(id, ptr);
		/* Read data */
		if constexpr (N > 0) {
			read_tuple<0, N>(t, ptr);
		}
		return true;		
	}

	/* Edge overload: read source and target ids first */
	template<typename T>
	bool read_tuple(size_t& source_id, size_t& target_id, T& t) {
		/* Size of tuple, can be zero */
		constexpr size_t N = std::tuple_size_v<T>;
		char* ptr = fgets(buffer, (int)bufferLength, file);
		if (ptr == nullptr) return false;
		/* Read source ID */
		ptr = Format< size_t >::read(source_id, ptr);
		ptr = Format< size_t >::read(target_id, ptr);
		/* Read data */
		if constexpr (N > 0) {
			read_tuple<0, N>(t, ptr);
		}
		return true;
	}
};