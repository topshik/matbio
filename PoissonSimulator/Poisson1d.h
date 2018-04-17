#pragma once
#include "stdafx.h"
#include "SplineBuilding.h"

using namespace std;
using namespace alglib;
using namespace spline_building;

#ifndef POISSON_1D_H
#define POISSON_1D_H

enum class BoundaryConditions_1d
{
	PERIODIC,
	KILLING,
	REFLECTIVE
};

struct Cell_1d
{

	vector<double> coords_x;
	vector<double> death_rates;

	Cell_1d() {}
};

template <typename RNG>
struct Grid_1d
{
	vector<Cell_1d> cells;
	vector<double> cell_death_rates;
	vector<int> cell_population;
	double Lx;
	int Nx;
	double b, d, dd;
	RNG rng;
	double initial_density;

	int total_population;
	double total_death_rate;

	double time;
	int event_count;

	double cutoff_radius_death;
	int cull_x;

	string death_kernel_string;
	spline1dinterpolant death_kernel_spline;

	string birth_kernel_string;
	spline1dinterpolant birth_reverse_cdf_spline;

	tuple<double, int> last_event;

	Cell_1d& cell_at(int i)
	{
		return cells[i];
	}

	double & cell_death_rate_at(int i)
	{
		return cell_death_rates[i];
	}

	int& cell_population_at(int i)
	{
		return cell_population[i];
	}


	void Initialize_death_rates()
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

	
	void kill_random()
	{

		int cell_death_index = boost::random::discrete_distribution<>(cell_death_rates)(rng);
		int in_cell_death_index = boost::random::discrete_distribution<>(cells[cell_death_index].death_rates)(rng);

		Cell_1d & death_cell = cells[cell_death_index];

		last_event = make_tuple<>(death_cell.coords_x[in_cell_death_index], -1);

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

	
	void spawn_random()
	{
		int cell_index = boost::random::discrete_distribution<>(cell_population)(rng);

		int event_index = boost::random::uniform_smallint<>(0, cell_population[cell_index] - 1)(rng);

		Cell_1d & parent_cell = cells[cell_index];

		double x_coord_new = parent_cell.coords_x[event_index] +
			spline1dcalc(birth_reverse_cdf_spline, boost::random::uniform_01<>()(rng))*(boost::random::bernoulli_distribution<>(0.5)(rng) * 2 - 1);

		last_event = make_tuple<>(x_coord_new, 1);

		if (x_coord_new<0 || x_coord_new>Lx)
		{
			last_event = make_tuple<>(x_coord_new, 0);
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

					double distance = abs(cell_at(i).coords_x[k] - x_coord_new);
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

	
	void make_event()
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

	
	void save_trajectory(ostream& output, double max_time, string type="") {
		for (auto cell : cells) {
			for (auto speciment_x : cell.coords_x) {
				output << fixed << setprecision(15) << "0" << "," << speciment_x << "," << "1" << "," << type << endl;
			}
		}
		while (time < max_time)
		{
			make_event();
			output << fixed << setprecision(15) << time << "," << get<0>(last_event) << "," << get<1>(last_event) << "," << type << endl;
		}
	}

	void save_trajectory(ostream& output, int max_events, string type = "") {
		for (auto cell : cells) {
			for (auto speciment_x : cell.coords_x) {
				output << fixed << setprecision(15) << "0" << "," << speciment_x << "," << "1" << "," << type << endl;
			}
		}
		for (int i = 0; i < max_events; i++) {
			make_event();
			output << fixed << setprecision(15) << time << "," << get<0>(last_event) << "," << get<1>(last_event) << "," << type << endl;
		}
	}

	void save_result(ostream& output, double max_time, string type = "") {
		while (time < max_time)
		{
			make_event();
		}

		for (auto cell : cells) {
			for (auto speciment_x : cell.coords_x) {
				output << fixed << setprecision(15) << speciment_x << "," << type << endl;
			}
		}

	}

	void save_result(ostream& output, int max_events, string type = "") {
		for (int i = 0; i < max_events; i++) {
			make_event();
		}

		for (auto cell : cells) {
			for (auto speciment_x : cell.coords_x) {
				output << fixed << setprecision(15) << speciment_x << "," << type << endl;
			}
		}
	}
	template <typename T>
	void save_slices(ostream& output, T params) {

		int discretization = get<1>(params);
		double stop_time= get<2>(params);
		double delta_time= get<3>(params);

		auto discretize = [&](int discretization) {
			vector<int> pops(discretization);
			for (auto cell : cells) {
				for (auto speciment_x : cell.coords_x) {
					int index = (int)floor(speciment_x / (Lx / discretization));
					if (index == discretization) index--;
					pops[index]++;
				}
			}
			return pops;
		};
		
		output << "#";
		boost::fusion::for_each(params, [&](auto item) {
			output << item << ",";
		});
		
		output << endl;

		output << "0";
		for (int pop : discretize(discretization)) {
			output << setprecision(15)<< "," << pop;
		}
		output << endl;

		int epochs = (int)floor(stop_time / delta_time);

		for (int i = 1; i <= epochs; i++) {
			while (time < delta_time*i) {
				make_event();
			}
			output << delta_time*i;
			for (int pop : discretize(discretization)) {
				output << "," << pop;
			}
			output << endl;
		}
	}

	Grid_1d(double Lx, int Nx, double b, double d, double dd, uint32_t Seed, double initial_density,
		string death_kernel, double death_cutoff, int death_spline_nodes,
		string birth_kernel, double birth_cutoff, int birth_spline_nodes) :
		cells(Nx), cell_death_rates(Nx), cell_population(Nx),
		Lx(Lx), Nx(Nx), b(b), d(d), dd(dd), rng(Seed),
		death_kernel_string(death_kernel),
		birth_kernel_string(birth_kernel),
		initial_density(initial_density), time(0), event_count(0), last_event(make_tuple<double, int>(0.0, 0))
	{
		cull_x = max(static_cast<int>(ceil(death_cutoff / (Lx / Nx))), 3);

		death_kernel_spline = Build_death_splines(death_kernel, death_cutoff, death_spline_nodes);
		birth_reverse_cdf_spline = Build_birth_splines(birth_kernel, birth_cutoff, birth_spline_nodes);

		Initialize_death_rates();
		total_death_rate = accumulate(cell_death_rates.begin(), cell_death_rates.end(), 0.0);
	}

};

#endif