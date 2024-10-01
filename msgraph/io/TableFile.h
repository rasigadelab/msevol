#pragma once
#pragma warning(disable:4996)

#include <cstdio>
#include <cstring>

class TableFile {
protected:
	/* Separator char */
	static constexpr char separator[] = "\t";
	/* File name */
	char* fname;
	/* File stream */
	FILE* file;

	inline void open(const char* fname_in, const char* access) {
		file = fopen(fname, access);
		/* HACK */
		if (file == nullptr) throw "Could not open file";
	}

	TableFile(const char* fname_in, const char* access) {
		fname = new char[strlen(fname_in) + 1];
		strcpy(fname, fname_in);
		open(fname, access);
	}

	TableFile(const char* path, const char* fname_in, const char* access) {
		constexpr const char* fnameFormat = "./%s/%s.csv";
		const size_t len =
			strlen(path) + strlen(fname_in) + strlen(fnameFormat);
		fname = new char[len + 1];
		sprintf(fname, fnameFormat, path, fname_in);
		open(fname, access);
	}

	~TableFile() {
		delete fname;
		if (file) fclose(file);
	}
};