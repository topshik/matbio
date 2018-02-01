#include "stdafx.h"

#include "Poisson1d.hpp"
#include "Poisson1d.cpp"

using namespace std;
using namespace boost;
using namespace chrono;
using namespace random;

/*
To run simulation, specify all parameters below and define CURRENT_DEATH_KERNEL(radius) and CULL_DEATH_INTERACTION_RADIUS
for death interaction function and define CURRENT_MOVE_RANDOM_VARIABLE_BY_ONE_AXIS(RNG) birth move random variable that
takes PRNG as a parameter.
*/
int main()
{
	double x_length_default = 1.0;
	int x_partition_default = 100;
	double birth_rate_default = 1e-4;
	double death_rate_default = 0;
	double comp_death_rate_default = 1e-4;
	uint32_t seed_default = 12345;
	double initial_density_default = 10.0;
	string death_kernel_default = "(1 / (sqrt(2*pi) * (2e-3)) * exp((-r^2) / (2 * (2e-3) * (2e-3))))";
	double death_interaction_cutoff_default = 6 * 2e-3;
	int death_spline_nodes_default = 1000;
	string birth_kernel_default = "(1 / (sqrt(2*pi) * (1e-2)) * exp((-r^2) / (2 * (1e-2) * (1e-2))))";
	double birth_interaction_cutoff_default = 6 * 1e-2;
	int birth_spline_death_nodes_default = 1000;

	double max_time_default = 300000;
	double time_step_default = 100;
	

	while (true) {
		stringstream file_name_stream;
		file_name_stream << "Simulation_" << duration_cast<seconds>(system_clock::now().time_since_epoch()).count() << ".txt";
		string file_name_default;
		file_name_stream >> file_name_default;
		string line;
		{
			cout << "Area length, press enter to keep default " << x_length_default << endl;
			getline(cin, line);
			if (!line.empty()) {
				istringstream(line) >> x_length_default;
				cout << "Area length is set to " << x_length_default << endl;
			}

			cout << "Partition count, press enter to keep default " << x_partition_default << endl;
			getline(cin, line);
			if (!line.empty()) {
				istringstream(line) >> x_partition_default;
				cout << "Partition count is set to " << x_partition_default << endl;
			}

			cout << "Birth rate, press enter to keep default " << birth_rate_default << endl;
			getline(cin, line);
			if (!line.empty()) {
				istringstream(line) >> birth_rate_default;
				cout << "Birth rate is set to " << birth_rate_default << endl;
			}

			cout << "Death rate, press enter to keep default " << death_rate_default << endl;
			getline(cin, line);
			if (!line.empty()) {
				istringstream(line) >> death_rate_default;
				cout << "Death rate is set to " << death_rate_default << endl;
			}

			cout << "Competitive death rate, press enter to keep default " << comp_death_rate_default << endl;
			getline(cin, line);
			if (!line.empty()) {
				istringstream(line) >> comp_death_rate_default;
				cout << "Competitive death rate is set to " << comp_death_rate_default << endl;
			}

			cout << "RNG seed, press enter to keep default " << seed_default << endl;
			getline(cin, line);
			if (!line.empty()) {
				istringstream(line) >> seed_default;
				cout << "RNG seed is set to " << seed_default << endl;
			}

			cout << "Initial density, press enter to keep default " << initial_density_default << endl;
			getline(cin, line);
			if (!line.empty()) {
				istringstream(line) >> initial_density_default;
				cout << "Initial density is set to " << initial_density_default << endl;
			}

			cout << "Death kernel, press enter to keep default " << death_kernel_default << endl;
			getline(cin, line);
			if (!line.empty()) {
				istringstream(line) >> death_kernel_default;
				cout << "Death kernel is set to " << death_kernel_default << endl;
			}

			cout << "Death interaction cutoff, press enter to keep default " << death_interaction_cutoff_default << endl;
			getline(cin, line);
			if (!line.empty()) {
				istringstream(line) >> death_interaction_cutoff_default;
				cout << "Death interaction cutoff is set to " << death_interaction_cutoff_default << endl;
			}

			cout << "Death kernel spline nodes, press enter to keep default " << death_spline_nodes_default << endl;
			getline(cin, line);
			if (!line.empty()) {
				istringstream(line) >> death_spline_nodes_default;
				cout << "Death kernel spline nodes count is set to " << death_spline_nodes_default << endl;
			}

			cout << "Birth kernel, press enter to keep default " << birth_kernel_default << endl;
			getline(cin, line);
			if (!line.empty()) {
				istringstream(line) >> birth_kernel_default;
				cout << "Birth kernel is set to " << birth_kernel_default << endl;
			}

			cout << "Birth interaction cutoff, press enter to keep default " << birth_interaction_cutoff_default << endl;
			getline(cin, line);
			if (!line.empty()) {
				istringstream(line) >> birth_interaction_cutoff_default;
				cout << "Birth interaction cutoff is set to " << birth_interaction_cutoff_default << endl;
			}

			cout << "Birth kernel spline nodes, press enter to keep default " << birth_spline_death_nodes_default << endl;
			getline(cin, line);
			if (!line.empty()) {
				istringstream(line) >> birth_spline_death_nodes_default;
				cout << "Birth kernel spline nodes count is set to " << birth_spline_death_nodes_default << endl;
			}

			cout << "Simulation time, press enter to keep default " << max_time_default << endl;
			getline(cin, line);
			if (!line.empty()) {
				istringstream(line) >> max_time_default;
				cout << "Simulation time is set to " << max_time_default << endl;
			}

			cout << "Data saving time step, press enter to keep default " << time_step_default << endl;
			getline(cin, line);
			if (!line.empty()) {
				istringstream(line) >> time_step_default;
				cout << "Simulation time is set to " << time_step_default << endl;
			}

			cout << "File name, press enter to keep default " << file_name_default << endl;
			getline(cin, line);
			if (!line.empty()) {
				istringstream(line) >> file_name_default;
				cout << "File name is set to " << file_name_default << endl;
			}
		}


		Grid_1d<lagged_fibonacci607> grid(
			x_length_default,					// Length by x
			x_partition_default,				// Partitions by x
			birth_rate_default,					// Birth rate
			death_rate_default,					// Death rate
			comp_death_rate_default,			// Competitive death rate
			seed_default,						// RNG seed
			initial_density_default,			// Initial density
			death_kernel_default,				// Death kernel
			death_interaction_cutoff_default,	// Death interaction cutoff
			death_spline_nodes_default,			// Death kernel spline nodes
			birth_kernel_default,				// Birth kernel
			birth_interaction_cutoff_default,	// Birth interaction cutoff
			birth_spline_death_nodes_default	// Birth kernel spline nodes
		);

		double i = 0;
		ofstream out(file_name_default);
		out << "Time,Population,Events" << endl;
		while (grid.time < max_time_default)
		{
			if (grid.time > i) {
				i += time_step_default;
				out << grid.time << "," << grid.total_population << "," << grid.event_count << endl;
				cout << grid.time << "," << grid.total_population << "," << grid.event_count << endl;
			}
			if (grid.total_population == 0) {
				cout << "Everything died";
				break;
			}
			grid.make_event();
		}
	}
	return 0;
}
