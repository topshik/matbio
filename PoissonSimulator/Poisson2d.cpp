#include "stdafx.h"
#include "Poisson2d.h"

#ifndef STARTING_CONDITIONS_2D
#define SIGMA_DEATH_2D 2e-3
#define SIGMA_MOVE_2D 1e-2

#define CULL_DEATH_INTERACTION_RADIUS_2D (3*SIGMA_DEATH_2D)
#define GAUSSIAN_KERNEL_2D(sigma,r2) (1 / (sqrt(2*M_PI) * (sigma)) * exp((-r2) / (2 * (sigma) * (sigma))))

#define CURRENT_DEATH_KERNEL_2D(r2) GAUSSIAN_KERNEL_2D(SIGMA_DEATH_2D,(r2))
#define CURRENT_MOVE_RANDOM_VARIABLE_BY_ONE_AXIS_2D(RNG) boost::random::normal_distribution<double>(0, SIGMA_MOVE_2D)(RNG)
#endif

using namespace std;
using namespace boost;
using namespace chrono;
using namespace random;

template <typename RNG>
void Grid_2d<RNG>::Initialize()
{
	//Spawn all speciments
	total_population = static_cast<int>(ceil(Lx*Ly*initial_density)); //initial population at t=0
	{
		double x_coord, y_coord;
		int i, j;
		for (int k = 0; k <total_population; k++)
		{
			x_coord = boost::uniform_real<>(0, Lx)(rng);
			y_coord = boost::uniform_real<>(0, Ly)(rng);

			i = static_cast<int>(floor(x_coord*Nx / Lx));
			j = static_cast<int>(floor(y_coord*Ny / Ly));

			if (i == Nx) i--;
			if (j == Ny) j--;

			cell_at(i, j).coords_x.push_back(x_coord);
			cell_at(i, j).coords_y.push_back(y_coord);

			cell_at(i, j).death_rates.push_back(d);
			cell_death_rate_at(i, j) += d;
			total_death_rate += d;

			cell_population_at(i, j)++;
		}
	}

	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			for (int k = 0; k < cell_population_at(i, j); k++)
			{
				for (int n = max(0, i - cull_x); n<min(Nx, i + cull_x + 1); n++)
				{
					for (int m = max(0, j - cull_y); m<min(Ny, j + cull_y + 1); m++)
					{
						for (int p = 0; p < cell_population_at(n, m); p++)
						{
							if (i == n && j == m && k == p) continue; // same speciment
							
							//Distance between k-th speciment in (i,j) cell and p-th speciment in (n,m) cell
							double delta_x = cell_at(i, j).coords_x[k] - cell_at(n, m).coords_x[p];
							double delta_y = cell_at(i, j).coords_y[k] - cell_at(n, m).coords_y[p];

							double interaction = dd*CURRENT_DEATH_KERNEL_2D(delta_x*delta_x + delta_y*delta_y);

							cell_at(i, j).death_rates[k] += interaction;
							cell_death_rate_at(i, j) += interaction;
							total_death_rate += interaction;
						}
					}
				}
			}
		}
	}
}

template <typename RNG>
void Grid_2d<RNG>::kill_random()
{

	int cell_death_index = boost::random::discrete_distribution<>(cell_death_rates)(rng);
	int in_cell_death_index = boost::random::discrete_distribution<>(cells[cell_death_index].death_rates)(rng);

	Cell_2d & death_cell = cells[cell_death_index];

	int cell_death_x = cell_death_index / Nx;
	int cell_death_y = cell_death_index % Nx;

	for (int i = max(0, cell_death_x - cull_x); i<min(Nx, cell_death_x + cull_x + 1); i++)
	{
		for (int j = max(0, cell_death_y - cull_y); j<min(Ny, cell_death_y + cull_y + 1); j++)
		{
			for (int k = 0; k < cell_population_at(i, j); k++)
			{
				if (i == cell_death_x && j == cell_death_y && k == in_cell_death_index) continue;

				double delta_x = cell_at(i, j).coords_x[k] - death_cell.coords_x[in_cell_death_index];
				double delta_y = cell_at(i, j).coords_y[k] - death_cell.coords_y[in_cell_death_index];

				double interaction = dd*CURRENT_DEATH_KERNEL_2D(delta_x*delta_x + delta_y*delta_y);

				cell_at(i, j).death_rates[k] -= interaction;
				//ignore dying speciment death rates since it is to be deleted

				cell_death_rate_at(i, j) -= interaction;
				cell_death_rate_at(cell_death_x, cell_death_y) -= interaction;

				total_death_rate -= 2 * interaction;
			}
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
	death_cell.coords_y[in_cell_death_index] = death_cell.coords_y[death_cell.coords_y.size() - 1];

	death_cell.death_rates.erase(death_cell.death_rates.end() - 1);
	death_cell.coords_x.erase(death_cell.coords_x.end() - 1);
	death_cell.coords_y.erase(death_cell.coords_y.end() - 1);
}

template <typename RNG>
void Grid_2d<RNG>::spawn_random()
{
	int cell_index = boost::random::discrete_distribution<>(cell_population)(rng);

	int event_index = boost::random::uniform_smallint<>(0, cell_population[cell_index] - 1)(rng);

	Cell_2d & parent_cell = cells[cell_index];

	double x_coord_new = parent_cell.coords_x[event_index] + CURRENT_MOVE_RANDOM_VARIABLE_BY_ONE_AXIS_2D(rng);
	double y_coord_new = parent_cell.coords_y[event_index] + CURRENT_MOVE_RANDOM_VARIABLE_BY_ONE_AXIS_2D(rng);

	if (x_coord_new<0 || x_coord_new>Lx || y_coord_new<0 || y_coord_new>Ly)
	{
		//Speciment failed to spawn and died outside area boundaries
	}
	else
	{

		int new_i = static_cast<int>(floor(x_coord_new*Nx / Lx));
		int new_j = static_cast<int>(floor(y_coord_new*Ny / Ly));

		if (new_i == Nx) new_i--;
		if (new_j == Ny) new_j--;

		//New speciment is added to the end of vector

		cell_at(new_i, new_j).coords_x.push_back(x_coord_new);
		cell_at(new_i, new_j).coords_y.push_back(y_coord_new);
		cell_at(new_i, new_j).death_rates.push_back(d);

		cell_death_rate_at(new_i, new_j) += d;
		total_death_rate += d;

		cell_population_at(new_i, new_j)++;
		total_population++;

		for (int i = max(0, new_i - cull_x); i<min(Nx, new_i + cull_x + 1); i++)
		{
			for (int j = max(0, new_j - cull_y); j<min(Ny, new_j + cull_y + 1); j++)
			{
				for (int k = 0; k < cell_population_at(i, j); k++)
				{
					if (i == new_i && j == new_j && k == cell_population_at(new_i, new_j) - 1) continue;

					double delta_x = cell_at(i, j).coords_x[k] - x_coord_new;
					double delta_y = cell_at(i, j).coords_y[k] - y_coord_new;

					double interaction = dd*CURRENT_DEATH_KERNEL_2D(delta_x*delta_x + delta_y*delta_y);

					cell_at(i, j).death_rates[k] += interaction;
					cell_at(new_i, new_j).death_rates[cell_population_at(new_i, new_j) - 1] += interaction;

					cell_death_rate_at(i, j) += interaction;
					cell_death_rate_at(new_i, new_j) += interaction;

					total_death_rate += 2 * interaction;
				}
			}
		}
	}
}

template <typename RNG>
void Grid_2d<RNG>::make_event()
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

template<typename RNG>
Grid_2d<RNG>::Grid_2d(double Lx, double Ly, int Nx, int Ny, double b, double d, double dd, uint32_t Seed, double initial_density):
	cells(Nx*Ny), cell_death_rates(Nx*Ny), cell_population(Nx*Ny),
	Lx(Lx), Ly(Ly), Nx(Nx), Ny(Ny), b(b), d(d), dd(dd), rng(Seed),
	initial_density(initial_density), time(0), event_count(0)
{
	cull_x = static_cast<int>(ceil(CULL_DEATH_INTERACTION_RADIUS_2D / (Lx / Nx)));
	cull_y = static_cast<int>(ceil(CULL_DEATH_INTERACTION_RADIUS_2D / (Ly / Ny)));

	Initialize();
	total_death_rate = accumulate(cell_death_rates.begin(), cell_death_rates.end(), 0.0);
}