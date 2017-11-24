#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <random>
#include <vector>

const double birth_rate = 0.4;
const double death_rate = 0.2;
double birth_variance = 0.1;
double death_variance = 0.1;
int wall = 0;

class Cell {
private:
    long  population, column;
    double x;
public:
    Cell() {};

    Cell(int Population, double X)
    : population(Population), x(X) {};

    long get_population() const {
        return population;
    }

    double get_coordinates() const {  // FIXME: for many dimensions need a vector
        return x;
    }

    long get_indices() const {  // FIXME: for many dimensions need a vector
        return column;
    }

    void set_population(long Population) {
        population = Population;
    }

    void set_coordinates(double X) {  // FIXME: for many dimensions need a vector
        x = X;
    }

    void set_indices(long Column) {  // FIXME: for many dimensions need a vector
        column = Column;
    }

    void print(std::ostream &out) const {
        out << std::setw(3) << population << std::endl;
    }

    friend std::ostream &operator<<(std::ostream &out, const Cell &cell) {
        cell.print(out);
        return out;
    }
};

class Grid {
private:
    long population, discretization;
    double width;  // FIXME: for many dimensions need a vector
    double cell_width;  // FIXME: for many dimensions need a vector
    std::vector<Cell> cells;
public:
    Grid() {};

    Grid(long Init_population, long Discretization, double Width)
    :population(Init_population), discretization(Discretization), width(Width) {
        cell_width = width / (double)discretization;
        cells = std::vector<Cell>(discretization);
        for (long i = 0; i != cells.size(); ++i) {
            cells[i].set_coordinates(width / discretization * i);
            cells[i].set_indices(i);
            cells[i].set_population(0);
        }
        for (long i = 0; i != population; ++i) {
            long index = (((long)rand() << 32) + rand()) % discretization;  // Слабоумие и слабоумие
            cells[index].set_population(cells[index].get_population() + 1);
        }
    }

    long get_population() const {
        return population;
    }

    void set_population(long new_population) {
        population = new_population;
    }

    long get_discretization() const {  // FIXME: for many dimensions need a vector
        return discretization;
    }

    double get_size() const {  // FIXME: for many dimensions need a vector
        return width;
    }

    double get_cell_size() const {  // FIXME: for many dimensions need a vector
        return cell_width;
    }

    std::vector<Cell> get_cells() const {  // FIXME: for many dimensions need a vector
        return cells;
    }

    Cell & operator[] (int x) {  // setter
        return cells[x];
    }

    Cell operator[] (int x) const {  // getter
        return cells[x];
    }
};

double distance(const Cell &from, const Cell &to) {  // FIXME: for many dimensions need a vector
    return std::sqrt(std::pow((from.get_coordinates() - to.get_coordinates()), 2));
}

// counting birth and death kernels (currently Gaussian)
double kernel(const Cell &cell1, const Cell &cell2, int type) {  // 0 - birth; 1 - death; FIXME: other kernels, many dim case
    if (!type) {
        return birth_rate / (sqrt(2 * M_PI * pow(birth_variance, 2))) * std::exp(-std::pow(distance(cell1, cell2), 2) / (2 * birth_variance));
    } else {
        return death_rate / (sqrt(2 * M_PI * pow(death_variance, 2))) * std::exp(-std::pow(distance(cell1, cell2), 2) / (2 * death_variance));
    }
}

// analysing which cells are approproate to be counted
std::vector<long> neighbour_birth_influence(const Grid &grid, const Cell &cell, int type) {  // 0 - birth; 1 - death
    double max_distance = 3 * birth_variance;  // 3 sigma rule
    if (type) max_distance = 3 * death_variance;  // 3 sigma rule
    long border_x = ceil(max_distance / grid.get_cell_size());
    std::vector<long> result;
    if (wall == 0) {  // killing border
        long left_border  = ((long)cell.get_indices() - border_x <= 0 ? 0 : cell.get_indices() - border_x);
        long right_border = ((long)cell.get_indices() + border_x >= grid.get_discretization() ? grid.get_discretization() : cell.get_indices() + border_x);
        for (long i = left_border; i != right_border; ++i) {
            result.push_back(i);
        }
    }
    return result;
}

std::vector<std::vector<long>> create_neighbours_vector(const Grid &grid, int type) {
    std::vector<std::vector<long>> result(grid.get_discretization());
    for (long i = 0; i != grid.get_discretization(); ++i) {
        result[i] = neighbour_birth_influence(grid, grid[i], type);
    }
    return result;
}

// life cycle of the grid
void iteration(Grid & grid, std::vector<std::vector<long>> &neighbours_birth, std::vector<std::vector<long>> &neighbours_death) {
    std::vector<double> nobirth_matrix(grid.get_discretization(), 1);
    std::vector<double> nodeath_matrix(grid.get_discretization(), 1);
    std::vector<long> cur_neighbours;
    for (long i = 0; i != grid.get_discretization(); ++i) {
        if (grid[i].get_population()) {
            cur_neighbours = neighbours_birth[i];
            for (int j : cur_neighbours) {
                double nobirth_prob = std::pow((1 - kernel(grid[i], grid[j], 0) * grid.get_cell_size()), grid[i].get_population());  // * cell_size instead of integration
                nobirth_matrix[grid[i].get_indices()] *= nobirth_prob;
            }
        }
    }
    for (long i = 0; i != grid.get_discretization(); ++i) {
        if (grid[i].get_population()) {
            cur_neighbours = neighbours_death[i];
            for (int j : cur_neighbours) {
                double nodeath_prob = std::pow((1 - kernel(grid[i], grid[j], 1) * grid.get_cell_size()), grid[i].get_population());  // * cell_size instead of integration
                nodeath_matrix[grid[j].get_indices()] *= nodeath_prob;
            }
        }
    }
    for (long i = 0; i != grid.get_discretization(); ++i) {
        if (grid[i].get_population() && (float)rand() / RAND_MAX < nodeath_matrix[i]) {
            grid[i].set_population(grid[i].get_population() - 1);
            grid.set_population(grid.get_population() - 1);
        }
        if ((float)rand() / RAND_MAX < nobirth_matrix[i]) {
            grid[i].set_population(grid[i].get_population() + 1);
            grid.set_population(grid.get_population() + 1);
        }
    }
}

int main(int argc, char ** argv) {
    srand(time(NULL));
    Grid grid(15000, 20000, 1);
    std::vector<std::vector<long>> neighbours_birth = create_neighbours_vector(grid, 0);
    std::vector<std::vector<long>> neighbours_death = create_neighbours_vector(grid, 1);
    for (int i = 0; i != 300; ++i) {
        std::cout << i << " " << grid.get_population();
        // for (int j = 0; j < grid.get_discretization(); ++j) {
        //     std::cout << " " << grid[j].get_population();
        // }
        std::cout << std::endl;
        iteration(grid, neighbours_birth, neighbours_death);
    }
}
