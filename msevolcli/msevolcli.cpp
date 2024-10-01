#include <iostream>
#include <cstdlib>
#include <filesystem>
#include <string>

#include "../msmodel/event_loop.h"
#include "../msmodel/universe_io.h"

/* Use four arguments: no. of steps between writes, no. of writes, input directory, output directory */

//using namespace std; 

/*
Expected arguments: none
*/
void write_reference_graph(const std::string outputDir)
{
	cout << "Writing reference graph to directory " << outputDir << std::endl;

	Universe u;
	std::filesystem::path path_out(outputDir);
	std::filesystem::create_directory(path_out);
	UniverseIO::writeGraph(u, outputDir.c_str(), { 0 }, false);
}

/*
Expected command arguments
SCHEMA refdir
RUN 10 10 test_in test_out 1 2
*/
int main(int argc, char* argv[])
{
	if (argc == 1) {
		std::cout <<
			"MSEvol command line interface v0.02\n" <<
			"Usage <COMMAND> <ARGUMENTS...>\n" <<
			"  SCHEMA <OUTPUT DIR> to write empty set of CSV files\n" <<
			"  RUN <STEPS BETWEEN WRITES> <WRITES> <INPUT DIR> <OUTPUT DIR> <SEED1> <SEED2> to run simulation\n";

		return EXIT_SUCCESS;
	}
	
std::string command(argv[1]);
	std::cout << command << std::endl;

	if (command == "SCHEMA") {
		write_reference_graph(argv[2]);
	}
	
	else if (command == "RUN") {

		size_t n_steps = atoll(argv[2]);
		size_t n_writes = atoll(argv[3]);

		std::string inputDir(argv[4]);
		std::string outputDir(argv[5]);

		// Random seeds
		uint32_t s1 = atoi(argv[6]);
		uint32_t s2 = atoi(argv[7]);

		std::cout << "No. of steps between writes: " << n_steps << std::endl;
		std::cout << "No. of writes: " << n_writes << std::endl;
		std::cout << "Input directory: " << inputDir << std::endl;
		std::cout << "Output directory: " << outputDir << std::endl;
		std::cout << flush;

		/* Universe setup */
		Universe u;
		UniverseIO::readGraph(u, inputDir.c_str());

		/* Serialize results (create directory if required) */
		std::filesystem::path path_out(outputDir);
		std::filesystem::create_directory(path_out);
		UniverseIO::writeGraph(u, outputDir.c_str(), { 0 }, false);

		/* Run events */
		Sampler sp(s1, s2);
		EventLoop ev(u, sp);
		for (size_t i = n_steps; i < n_steps * (n_writes + 1); i += n_steps) {
			ev.run(n_steps);
			UniverseIO::writeGraph(u, outputDir.c_str(), { i }, true);
			auto s = u.size();
			std::cout << "Iter = " << i
				<< ", |V| = " << s.first
				<< " --- |E| = " << s.second << std::endl << flush;
		}
		std::cout << std::endl;

		return EXIT_SUCCESS;

	}
	else {
		std::cout << "Bad arguments, exiting." << std::endl;
		return EXIT_FAILURE;
	}
}
