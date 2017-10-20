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
int wall = 0;

class Cell {
private:
    unsigned long long  population, column;
    double x;
public:
    Cell() {};

    Cell(int Population, double X)
    : population(Population), x(X) {};

    unsigned long long get_population() const {
        return population;
    }

    double get_coordinates() const {  // FIXME: for many dimensions need a vector
        return x;
    }

    unsigned long long get_indices() const {  // FIXME: for many dimensions need a vector
        return column;
    }

    void set_population(unsigned long long Population) {
        population = Population;
    }

    void set_coordinates(double X) {  // FIXME: for many dimensions need a vector
        x = X;
    }

    void set_indices(unsigned long long Column) {  // FIXME: for many dimensions need a vector
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
    unsigned long long population, discretization;
    double width;  // FIXME: for many dimensions need a vector
    std::vector<Cell> cells;
    double cell_width;  // FIXME: for many dimensions need a vector
public:
    Grid() {};

    Grid(unsigned long long Init_population, unsigned long long Discretization, double Width)
    :population(Init_population), discretization(Discretization), width(Width) {
        cells = std::vector<Cell>(discretization);
        Cell init_cell;
        for (unsigned long long i = 0; i != population; ++i) {
            init_cell.set_population(1);
            init_cell.set_indices(rand() % discretization);             // FIXME: indices, columns, rows??
            init_cell.set_coordinates((width / discretization) * ((double)init_cell.get_indices() + 0.5));
            cells[init_cell.get_indices()] = init_cell;
        }
    }

    unsigned long long get_population() const {
        return population;
    }

    unsigned long long get_discretization() const {  // FIXME: for many dimensions need a vector
        return discretization;
    }

    double get_size() const {  // FIXME: for many dimensions need a vector
        return width;
    }

    double get_cell_size() const {  // FIXME: for many dimensions need a vector
        return cell_width;
    }

    std::vector<Cell> get_cells() const {  // FIXME: for many dimensions need a vector
        return cells;  // not sure about const to edit with cell.set_value()
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

std::vector<Cell> neighbour_birth_influence(const Grid &grid, const Cell &cell) {
    // max_distance ? need func for arbitrary Gaussian distribution expectation and variance
    double max_distance = 0.5;
    unsigned long long border_x = ceil(max_distance / grid.get_cell_size());  // check for types conversion
    std::vector<Cell> result;
    if (wall == 0) {
        for (unsigned long long i = (cell.get_indices() - border_x < 0 ? 0 : cell.get_indices() - border_x);
                                i != (cell.get_indices() + border_x >= grid.get_discretization() ? grid.get_discretization() : cell.get_indices() + border_x);
                                ++i) {
            result.push_back(grid.get_cells()[i]);
        }
    }
}

int main(int argc, char ** argv) {
    srand(time(NULL));
    Grid grid(10, 10, 1);
    for (auto cell : neighbour_birth_influence(grid, grid[0]))
        std::cout << cell;
}
