#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <random>
#include <vector>

const double birth_rate = 0.4;
const double death_rate = 0.2;

class Cell {
public:
    int population, row, column;
    double x, y;
    Cell() {};
    Cell(int Population, double X, double Y)
    : population(Population), x(X), y(Y) {};
};

class Grid {
public:
    double width, height;
    size_t n, total_population;         // FIXME: private and public fields
    std::vector<std::vector<Cell> > cells;
    Grid() {};
    Grid(size_t Width, size_t Height, size_t N, size_t Init_population)
    : width(Width), height(Height), n(N), total_population(Init_population) {
        cells = std::vector<std::vector<Cell> >(n, std::vector<Cell>(n));
        Cell init_cell;

        for (size_t i = 0; i != total_population; ++i) {
            init_cell.population = 1;
            init_cell.row = rand() % n;
            init_cell.column = rand() % n;
            init_cell.x = (width / n) * ((double)init_cell.column + 0.5);
            init_cell.y = (height / n) * ((double)init_cell.row + 0.5);
            cells[init_cell.row][init_cell.column] = init_cell;
        }
    };

    void print(std::ostream &out) const {
        for (auto i : cells) {
            for (auto j : i) {
                out << std::setw(3) << j.population << " ";
            }
            out << std::endl;
        }
    }

    friend std::ostream &operator<<(std::ostream &out, const Grid &grid) {
        grid.print(out);
        return out;
    }

    std::vector<Cell> & operator[] (int x) {
        // std::cout << "post" << std::endl;             // Debug line
        return cells[x];
    }
    std::vector<Cell> operator[] (int x) const {
        // std::cout << "get" << std::endl;              // Debug line
        return cells[x];
    }
};

double distance(const Cell &from, const Cell &to) {
    return std::sqrt(std::pow((from.x - to.x), 2) + std::pow((from.y - to.y), 2));
}

/* Analyzing which cells are appropriate to be counted */
double max_distance = 3;
std::vector<Cell> neighbour_birth_influence(const Grid &grid, const Cell &cell) {
    std::vector<Cell> result;
    int border_x = ceil(max_distance * grid.n / grid.width);  // Max distance from cell by x
    int border_y = ceil(max_distance * grid.n / grid.height);  // Max distance from cell by y
    for (size_t i = (cell.column - border_x < 0 ? 0 : cell.column - border_x); i <
                    (cell.column + border_x >= grid.n ? grid.n : cell.column + border_x); ++i) {
        for (size_t j = (cell.row - border_y < 0 ? 0 : cell.row - border_y); j <
                    (cell.row + border_y >= grid.n ? grid.n : cell.row + border_y); ++j) {
            if (distance(cell, grid.cells[i][j]) < max_distance) {
                result.push_back(grid.cells[i][j]);
            }
        }
    }
    return result;
};

void iteration(Grid &grid) {
    std::vector<std::vector<double> > nobirth_matrix(grid.n, std::vector<double>(grid.n, 1));
    /* Counting (no)birth probabilities matrix*/
    for (auto row : grid.cells) {
        for (auto cell : row) {
            if (cell.population > 0) {
                std::vector<Cell> neighbours = neighbour_birth_influence(grid, cell);
                for (auto neighbour : neighbours) {
                    double cell_influence = pow(std::exp(-std::pow(distance(neighbour, cell), 2) / 2), cell.population);
                    double birth_prob = birth_rate * cell_influence;
                    nobirth_matrix[neighbour.row][neighbour.column] *= (1 - birth_prob);
                }
            }
        }
    }

    /* Giving birth and killing*/
    for (size_t i = 0; i != nobirth_matrix.size(); ++i) {
        for (size_t j = 0; j != nobirth_matrix[0].size(); ++j) {
            if ((double)rand() / RAND_MAX > nobirth_matrix[i][j]) {
                ++grid.cells[i][j].population;
                ++grid.total_population;
            }

            int died = 0;
            for (size_t k = 0; k != grid.cells[i][j].population; ++k) {
                if ((double)rand() / RAND_MAX < death_rate) {
                    ++died;
                }
            }

            grid.cells[i][j].population -= died;
            grid.total_population -= died;
        }
    }
}

int main(int argc, char **argv) {

    std::string usage_string = "Usage: " + std::string(argv[0]) + " [size] [discretization] [iterations] [initial population]";
    if (argc != 5) {
        std::cerr << "Wrong number of arguments\n" << usage_string << std::endl;
        return 1;
    }
    char *endptr;
    long int size = strtol(argv[1], &endptr, 10);
    if (!*argv[1] || *endptr) {
        std::cerr << "Wrong size: " << argv[1] << std::endl << usage_string << std::endl;
        return 1;
    }
    long int discretization = strtol(argv[2], &endptr, 10);
    if (!*argv[2] || *endptr) {
        std::cerr << "Wrong discretization: " << argv[2] << std::endl << usage_string << std::endl;
        return 1;
    }
    long int iterations = strtol(argv[3], &endptr, 10);
    if (!*argv[3] || *endptr) {
        std::cerr << "Wrong number of iterations: " << argv[3] << std::endl << usage_string << std::endl;
        return 1;
    }
    long int init_population = strtol(argv[4], &endptr, 10);
    if (!*argv[4] || *endptr) {
        std::cerr << "Wrong initial population: " << argv[4] << std::endl;
        return 1;
    }

    srand(time(NULL));
    Grid grid(size, size, discretization, init_population);
    std::cout << 0 << "\t" << grid.total_population << std::endl;
    for (long int i = 1; i <= iterations; ++i) {
        iteration(grid);
        std::cout << i << "\t" << grid.total_population << std::endl;
    }
    return 0;
}
