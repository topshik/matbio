#include "stdafx.h"
#include "Poisson1d.hpp"


using namespace std;
using namespace boost;
using namespace chrono;
using namespace random;
using namespace alglib;

	template <typename RNG>
	void Grid_1d<RNG>::Build_death_splines(string death_kernel, double death_cutoff, int death_spline_nodes) {

	double death_r;

	exprtk::symbol_table<double>	death_kernel_symbol_table;
	exprtk::expression<double>		death_kernel_expression;
	exprtk::parser<double>			death_kernel_parser;

	death_kernel_symbol_table.add_variable("r", death_r);
	death_kernel_symbol_table.add_constants();
	death_kernel_expression.register_symbol_table(death_kernel_symbol_table);
	death_kernel_parser.compile(death_kernel_string, death_kernel_expression);

	real_1d_array x_1d_array;
	real_1d_array y_1d_array;

	x_1d_array.setlength(death_spline_nodes + 3);
	y_1d_array.setlength(death_spline_nodes + 3);

	for (int i = 0; i < death_spline_nodes + 3; i++) {
		x_1d_array[i] = death_cutoff*i / (death_spline_nodes + 3 - 1);
		death_r = x_1d_array[i];
		if (i < death_spline_nodes){
			y_1d_array[i] = death_kernel_expression.value();
		}
		else {
			y_1d_array[i] = 0;
		}
	}

	spline1dbuildmonotone(x_1d_array, y_1d_array, death_kernel_spline);
}

	template <typename RNG>
	void Grid_1d<RNG>::Build_birth_splines(string birth_kernel, double birth_cutoff, int birth_spline_nodes) {
		double birth_r;

		exprtk::symbol_table<double>	birth_kernel_symbol_table;
		exprtk::expression<double>		birth_kernel_expression;
		exprtk::parser<double>			birth_kernel_parser;

		birth_kernel_symbol_table.add_variable("r", birth_r);
		birth_kernel_symbol_table.add_constants();
		birth_kernel_expression.register_symbol_table(birth_kernel_symbol_table);
		birth_kernel_parser.compile(birth_kernel_string, birth_kernel_expression);

		real_1d_array x_1d_array;
		real_1d_array y_1d_array;

		x_1d_array.setlength(birth_spline_nodes+3);
		y_1d_array.setlength(birth_spline_nodes+3);

		for (int i = 0; i < birth_spline_nodes+3; i++) {
			x_1d_array[i] = birth_cutoff*i / (birth_spline_nodes +3 - 1);
			birth_r = x_1d_array[i];
			if (i < birth_spline_nodes) {
				y_1d_array[i] = birth_kernel_expression.value();
			}
			else {
				y_1d_array[i] = 0;
			}
		}
		spline1dbuildmonotone(x_1d_array, y_1d_array, birth_kernel_spline);

		real_1d_array x_quantile_1d_array;
		real_1d_array y_quantile_1d_array;

		x_quantile_1d_array.setlength(birth_spline_nodes);
		y_quantile_1d_array.setlength(birth_spline_nodes);

		double approx_const = spline1dintegrate(birth_kernel_spline, birth_cutoff);
		for (int i = 0; i < birth_spline_nodes; i++) {
			x_quantile_1d_array[i] = (double)i / (birth_spline_nodes - 1);
			y_quantile_1d_array[i] = 
				math::tools::newton_raphson_iterate(
					[=](double y) {return make_tuple(
						spline1dintegrate(birth_kernel_spline, y)/approx_const - x_quantile_1d_array[i],
						spline1dcalc(birth_kernel_spline,y) / approx_const); },
				1e-10, 0.0, birth_cutoff, numeric_limits<double>::digits);
		}
		spline1dbuildmonotone(x_quantile_1d_array, y_quantile_1d_array, birth_reverse_cdf_spline);
	}
	

	template <typename RNG>
	void Grid_1d<RNG>::Initialize_death_rates()
	{
		//Spawn all speciments
		total_population = static_cast<int>(ceil(Lx*initial_density)); //initial population at t=0
		{
			double x_coord;
			int i;
			for (int k = 0; k < total_population; k++)
			{
				x_coord = boost::uniform_real<>(0, Lx)(rng);

				i = static_cast<int>(floor(x_coord*Nx / Lx));

				if (i == Nx) i--;

				cell_at(i).coords_x.push_back(x_coord);

				cell_at(i).death_rates.push_back(d);
				cell_death_rate_at(i) += d;
				total_death_rate += d;

				cell_population_at(i)++;
			}
		}

		for (int i = 0; i < Nx; i++)
		{
			for (int k = 0; k < cell_population_at(i); k++)
			{
				for (int n = max(0, i - cull_x); n < min(Nx, i + cull_x + 1); n++)
				{
					for (int p = 0; p < cell_population_at(n); p++)
					{
						if (i == n && k == p) continue; // same speciment

														//Distance between k-th speciment in (i) cell and p-th speciment in (n) cell

						double distance = abs(cell_at(i).coords_x[k] - cell_at(n).coords_x[p]);

						double interaction = dd*spline1dcalc(death_kernel_spline, distance);

						cell_at(i).death_rates[k] += interaction;
						cell_death_rate_at(i) += interaction;
						total_death_rate += interaction;
					}
				}
			}
		}
	}

	template <typename RNG>
	void Grid_1d<RNG>::kill_random()
	{

		int cell_death_index = boost::random::discrete_distribution<>(cell_death_rates)(rng);
		int in_cell_death_index = boost::random::discrete_distribution<>(cells[cell_death_index].death_rates)(rng);

		Cell_1d & death_cell = cells[cell_death_index];

		int cell_death_x = cell_death_index;

		for (int i = max(0, cell_death_x - cull_x); i < min(Nx, cell_death_x + cull_x + 1); i++)
		{
			for (int k = 0; k < cell_population_at(i); k++)
			{
				if (i == cell_death_x && k == in_cell_death_index) continue;

				
				double distance = abs(cell_at(i).coords_x[k] - death_cell.coords_x[in_cell_death_index]);
				double interaction = dd*spline1dcalc(death_kernel_spline, distance);

				cell_at(i).death_rates[k] -= interaction;
				//ignore dying speciment death rates since it is to be deleted

				cell_death_rate_at(i) -= interaction;
				cell_death_rate_at(cell_death_x) -= interaction;

				total_death_rate -= 2 * interaction;
			}
		}
		//remove dead speciment
		cell_death_rates[cell_death_index] -= d;
		total_death_rate -= d;

		if (abs(cell_death_rates[cell_death_index])<1e-10)
		{
			cell_death_rates[cell_death_index] = 0;
		}

		cell_population[cell_death_index]--;
		total_population--;

		//swap dead and last
		death_cell.death_rates[in_cell_death_index] = death_cell.death_rates[death_cell.death_rates.size() - 1];
		death_cell.coords_x[in_cell_death_index] = death_cell.coords_x[death_cell.coords_x.size() - 1];

		death_cell.death_rates.erase(death_cell.death_rates.end() - 1);
		death_cell.coords_x.erase(death_cell.coords_x.end() - 1);
	}

	template <typename RNG>
	void Grid_1d<RNG>::spawn_random()
	{
		int cell_index = boost::random::discrete_distribution<>(cell_population)(rng);

		int event_index = boost::random::uniform_smallint<>(0, cell_population[cell_index] - 1)(rng);

		Cell_1d & parent_cell = cells[cell_index];

		double x_coord_new = parent_cell.coords_x[event_index] + 
			spline1dcalc(birth_reverse_cdf_spline, boost::random::uniform_01<>()(rng))*(boost::random::bernoulli_distribution<double>(0.5)(rng)*2-1);

		if (x_coord_new<0 || x_coord_new>Lx)
		{
			//Speciment failed to spawn and died outside area boundaries
		}
		else
		{

			int new_i = static_cast<int>(floor(x_coord_new*Nx / Lx));

			if (new_i == Nx) new_i--;

			//New speciment is added to the end of vector

			cell_at(new_i).coords_x.push_back(x_coord_new);
			cell_at(new_i).death_rates.push_back(d);

			cell_death_rate_at(new_i) += d;
			total_death_rate += d;

			cell_population_at(new_i)++;
			total_population++;

			for (int i = max(0, new_i - cull_x); i < min(Nx, new_i + cull_x + 1); i++)
			{
				for (int k = 0; k < cell_population_at(i); k++)
				{
					if (i == new_i && k == cell_population_at(new_i) - 1) continue;

					double distance =abs(cell_at(i).coords_x[k] - x_coord_new);
					double interaction = dd*spline1dcalc(death_kernel_spline, distance);

					cell_at(i).death_rates[k] += interaction;
					cell_at(new_i).death_rates[cell_population_at(new_i) - 1] += interaction;

					cell_death_rate_at(i) += interaction;
					cell_death_rate_at(new_i) += interaction;

					total_death_rate += 2 * interaction;
				}
			}
		}

	}

	template <typename RNG>
	void Grid_1d<RNG>::make_event()
	{
		event_count++;
		time += boost::random::exponential_distribution<>(total_population*b + total_death_rate)(rng);

		if (boost::random::bernoulli_distribution<>(total_population*b / (total_population*b + total_death_rate))(rng) == 0) //Someone died
		{
			kill_random();
		}
		else //Someone was born. 
		{
			spawn_random();
		}

	}

	template <typename RNG>
	Grid_1d<RNG>::Grid_1d(double Lx, int Nx, double b, double d, double dd, uint32_t Seed, double initial_density,
						  string death_kernel, double death_cutoff, int death_spline_nodes,
						  string birth_kernel, double birth_cutoff, int birth_spline_nodes) :
		cells(Nx), cell_death_rates(Nx), cell_population(Nx),
		Lx(Lx), Nx(Nx), b(b), d(d), dd(dd), rng(Seed),
		death_kernel_string(death_kernel), death_kernel_spline(),
		birth_kernel_string(birth_kernel), birth_kernel_spline(),birth_reverse_cdf_spline(),
		initial_density(initial_density), time(0), event_count(0)
	{
		cull_x = max(static_cast<int>(ceil(death_cutoff / (Lx / Nx))),3);

		Build_death_splines(death_kernel, death_cutoff, death_spline_nodes);
		Build_birth_splines(birth_kernel, birth_cutoff, birth_spline_nodes);
		Initialize_death_rates();
		total_death_rate = accumulate(cell_death_rates.begin(), cell_death_rates.end(), 0.0);
	}