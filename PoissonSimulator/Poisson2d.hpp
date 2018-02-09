#pragma once
#include "stdafx.h"
using namespace std;

#ifndef POISSON_2D_H
#define POISSON_2D_H

enum class BoundaryConditions_2d
{
	PERIODIC,
	KILLING,
	REFLECTIVE
};

struct Cell_2d
{

	vector<double> coords_x;
	vector<double> coords_y;
	vector<double> death_rates;

	Cell_2d() {}
};

template <typename RNG>
struct Grid_2d
{
	vector<Cell_2d> cells;
	vector<double> cell_death_rates;
	vector<int> cell_population;
	double Lx, Ly;
	int Nx, Ny;
	double b, d, dd;
	RNG rng;
	double initial_density;

	int total_population;
	double total_death_rate;

	double time;
	int event_count;

	int cull_x;
	int cull_y;

	Cell_2d& cell_at(int i, int j)
	{
		return cells[Nx*i + j];
	}

	double & cell_death_rate_at(int i, int j)
	{
		return cell_death_rates[Nx*i + j];
	}

	int& cell_population_at(int i, int j)
	{
		return cell_population[Nx*i + j];
	}

	void Initialize();
	void kill_random();
	void spawn_random();
	void make_event();

	Grid_2d(double Lx, double Ly, int Nx, int Ny, double b, double d, double dd, uint32_t Seed, double initial_density);
};

#endif
