#include "stdafx.h"

#include "Poisson1d.h"

using namespace std;
using namespace boost;
using namespace chrono;
using namespace random;
/*

*/
int main(int argc, char **argv)
{

	/*
	size discretization stop_time delta_time initial_population birth_rate death_rate birth_variance death_variance[seed]

	parameter example:
	1 10000 1000 10 500 6e-3 1e-5 2e-1 0.001 2
	*/

	double x_length = atof(argv[1]);
	int discretization = atoi(argv[2]);
	double stop_time = atof(argv[3]);
	double delta_time = atof(argv[4]);
	int initial_population = atoi(argv[5]);
	double birth_rate= atof(argv[6]);
	double death_rate= atof(argv[7]);
	double birth_var= atof(argv[8]);
	double death_var= atof(argv[9]);
	int seed= atoi(argv[10]);

	auto params= make_tuple(x_length, discretization, stop_time, delta_time, initial_population, birth_rate, death_rate, birth_var,death_var,seed);
	
	string death_kernel;
	string birth_kernel;

	stringstream str_str;
	str_str << "(1 / (sqrt(2*pi) * " << birth_var << ") * exp((-r^2) / (2 * " << birth_var << " * " << birth_var << ")))" << endl;
	getline(str_str, birth_kernel);

	str_str << "(1 / (sqrt(2*pi) * " << death_var << ") * exp((-r^2) / (2 * " << death_var << " * " << death_var << ")))" << endl;
	getline(str_str, death_kernel);

	Grid_1d<lagged_fibonacci607> grid(
		x_length,								// Length by x
		250,									// Partitions by x
		birth_rate,								// Birth rate
		0,										// Death rate
		death_rate,								// Competitive death rate
		seed,									// RNG seed
		(double)initial_population / x_length,	// Initial density
		death_kernel,							// Death kernel
		death_var*6,							// Death interaction cutoff
		1000,									// Death kernel spline nodes
		birth_kernel,							// Birth kernel
		birth_var * 6,							// Birth interaction cutoff
		1000									// Birth kernel spline nodes
	);

	grid.save_slices<decltype(params)>(cout, params);
	/*
	double x_length_default = 1.0;
	int x_partition_default = 100;
	double birth_rate_default = 1e-4;
	double death_rate_default = 0;
	double comp_death_rate_default = 1e-4;
	uint32_t seed_default = 12345;
	double initial_density_default = 10.0;
	string death_kernel_default = "(1 / (sqrt(2*pi) * (1e-2)) * exp((-r^2) / (2 * (1e-2) * (1e-2))))";
	double death_interaction_cutoff_default = 6 * 1e-2;
	int death_spline_nodes_default = 1000;
	string birth_kernel_default = "(1 / (sqrt(2*pi) * (2e-2)) * exp((-r^2) / (2 * (2e-2) * (2e-2))))";
	double birth_interaction_cutoff_default = 6 * 2e-2;
	int birth_spline_death_nodes_default = 1000;

	double max_time_default = 100;
	double time_step_default = 100;
	


	birth_rate_default = 6.158482;
	comp_death_rate_default = 1;
	birth_kernel_default = "(1 / (sqrt(2*pi) * 0.0112883789168469) * exp((-r^2) / (2 * 0.0112883789168469 * 0.0112883789168469)))";
	birth_interaction_cutoff_default = 6 * 0.0112883789168469;
	x_length_default = 25;

	vector<double> sigma_death = {	0.001,				0.00127427498570313,0.00162377673918872,0.00206913808111479,0.00263665089873036,
									0.00335981828628378,0.00428133239871939,0.00545559478116852,0.00695192796177561,0.00885866790410083,
									0.0112883789168469,	0.0143844988828766,	0.0183298071083244,	0.0233572146909012,	0.0297635144163132,
									0.0379269019073225,	0.0483293023857175,	0.0615848211066026,	0.0784759970351461,	0.1 };
	vector<double> death_cutoffs(20, 0);
	vector<string> death_kernels(20,"");

	vector<double> x_multipliers = {	1.14746328357921,1.42709906070841,1.77808737346116,2.21969384653463,2.79454291002368,
										3.54740592420392,4.66439949930841,6.39874063436738,9.58726933921289,17.1577302445848,
										32.107358756701,64.4361606200025,100.555185248586,169.192515667306,207.935477448126,
										219.130647312954,213.957775386685,219.296668742568,177.693727966807,162.776505265562 };
	for (auto i = 0; i < 20; i++) {
		death_cutoffs[i] = 6 * sigma_death[i];

		stringstream str_str;
		str_str << "(1 / (sqrt(2*pi) * " << sigma_death[i] << ") * exp((-r^2) / (2 * " << sigma_death[i] << " * " << sigma_death[i] << ")))" << endl;
		getline(str_str,death_kernels[i]);
	}


	for (int i = 0; i < 20; i++) {
		stringstream file_name_stream;
		file_name_stream << "Simulation_" << i << ".csv";
		string file_name;
		file_name_stream >> file_name;
		ofstream out(file_name);

		Grid_1d<lagged_fibonacci607> grid(
			x_length_default*x_multipliers[i],	// Length by x
			x_partition_default,				// Partitions by x
			birth_rate_default,					// Birth rate
			death_rate_default,					// Death rate
			comp_death_rate_default,			// Competitive death rate
			seed_default,						// RNG seed
			5000.0/(x_length_default*x_multipliers[i]),// Initial density
			death_kernels[i],					// Death kernel
			death_cutoffs[i],					// Death interaction cutoff
			death_spline_nodes_default,			// Death kernel spline nodes
			birth_kernel_default,				// Birth kernel
			birth_interaction_cutoff_default,	// Birth interaction cutoff
			birth_spline_death_nodes_default	// Birth kernel spline nodes
		);
		grid.save_trajectory(out, 500000*(i+1),to_string(i));
	}

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

		grid.save_trajectory(out, max_time_default*100);
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
	*/
	return 0;
}
