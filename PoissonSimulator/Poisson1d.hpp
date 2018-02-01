#pragma once
#include "stdafx.h"
#include "exprtk.hpp"
#include "interpolation.h"

using namespace std;
using namespace alglib;

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
	spline1dinterpolant birth_kernel_spline;
	spline1dinterpolant birth_reverse_cdf_spline;

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

	void Initialize_death_rates();

	void Build_death_splines(string death_kernel, double death_cutoff, int death_spline_nodes);
	void Build_birth_splines(string birth_kernel, double birth_cutoff, int birth_spline_nodes);

	void kill_random();
	void spawn_random();
	void make_event();

	Grid_1d(double Lx, int Nx, 
		double b, double d, double dd, 
		uint32_t Seed, double initial_density,
		string death_kernel, double death_cutoff, int death_spline_nodes,
		string birth_kernel, double birth_cutoff, int birth_spline_nodes
	);
};

#endif
