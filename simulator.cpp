#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <random>
#include <vector>

const double birth_rate = 0.4;
const double death_rate = 0;
int wall;

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
            init_cell.x = ((double)width / n) * ((double)init_cell.column + 0.5);
            init_cell.y = ((double)height / n) * ((double)init_cell.row + 0.5);
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
std::vector<Cell> neighbour_birth_influence(const Grid &grid, const Cell &cell) {  // 
    std::vector<Cell> result;
    int border_x = ceil(max_distance * grid.n / grid.width);  // Max distance from cell by x
    int border_y = ceil(max_distance * grid.n / grid.height);  // Max distance from cell by y
    if (wall == 0) {
        for (size_t i = (cell.row - border_y < 0 ? 0 : cell.row - border_y); i <
                        (cell.row + border_y >= (int)grid.n ? grid.n : cell.row + border_y); ++i) {
            for (size_t j = (cell.column - border_x < 0 ? 0 : cell.column - border_x); j <
                            (cell.column + border_x >= (int)grid.n ? grid.n : cell.column + border_x); ++j) {
                if (distance(cell, grid.cells[i][j]) < max_distance) {
                    result.push_back(grid.cells[i][j]);
                }
            }
        }
    } else if (wall == 1) {
        for (int i = cell.row - border_y; i < cell.row + border_y; ++i) {
            for (int j = cell.column - border_x; j < cell.column + border_x; ++j) {
                int row, col;
            
                if (cell.column - border_x < 0) {
                    col = -j;
                } else if (cell.column + border_x >= (int)grid.n) {
                    col = 2 * grid.n - 2 - j;
                } else {
                    col = j;
                }

                if (cell.row - border_y < 0) {
                    row = -i;
                } else if (cell.row + border_y >= (int)grid.n) {
                    row = 2 * grid.n - 2 - i    ;
                } else {
                    row = i;
                }

                if (col >= 0 && col < (int)grid.n && row >= 0 && row < (int)grid.n &&
                    distance(cell, grid.cells[row][col]) < max_distance) {
                    result.push_back(grid.cells[row][col]);
                }
            }
        }
    } else { 
        for (int i = cell.row - border_y; i < cell.row + border_y; ++i) {
            for (int j = cell.column - border_x; j < cell.column + border_x; ++j) {
                int row, col;

                if (cell.column - border_x < 0) {
                    col = 0;
                } else if (cell.column + border_x >= (int)grid.n) {
                    col = (int)grid.n - 1;
                } else {
                    col = j;
                }

                if (cell.row - border_y < 0) {
                    row = 0;
                } else if (cell.row + border_y >= (int)grid.n) {
                    row = (int)grid.n - 1;
                } else {
                    row = i;
                }

                if (distance(cell, grid.cells[row][col]) < max_distance) {
                    result.push_back(grid.cells[row][col]);
                }
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
                    double cell_influence = pow(std::exp(-std::pow(distance(neighbour, cell), 2) / 2), cell.population); // TODO: No pow
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
            for (int k = 0; k != grid.cells[i][j].population; ++k) {
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

    std::string usage_string = "Usage: " + std::string(argv[0]) + " [size] [discretization] [iterations] [initial population] [wall type]";
    if (argc != 6) {
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

    wall = strtol(argv[5], &endptr, 10);
    if (!*argv[5] || *endptr) {
        std::cerr << "Wrong wall type:" << argv[5] << "\n0 - killing, 1 - reflecting, 2 - reflecting to wall"<< std::endl << usage_string << std::endl;
        return 1;
    }

    std::ofstream population, moment1, init_density, density;

    srand(42);
    Grid grid(size, size, discretization, init_population);

    init_density.open("init_density.csv");
    for (auto row : grid.cells) {
        for (int i = 0; i < (int)row.size(); ++i) {
            init_density << row[i].population;
            if (i < (int)row.size() - 1) {
                init_density << ',';
            }
        }
        init_density << std::endl;
    }
    init_density.close();

    population.open("population.csv");
    moment1.open("moment1.csv");
    std::cout << 0 << "\t" << grid.total_population << std::endl;
    for (long int i = 1; i <= iterations; ++i) {
        iteration(grid);
        std::cout << i << "\t" << grid.total_population << std::endl;
        population << grid.total_population;
        moment1 << (double)grid.total_population / (grid.width * grid.height);
        if (i != iterations) {
            population << ',';
            moment1 << ',';
        }
        if (i % 10 == 1) {
            density.open("density" + std::to_string(i) + ".csv");
            for (auto row : grid.cells) {
                for (int j = 0; j < (int)row.size(); ++j) {
                    density << row[j].population;
                    if (j < (int)row.size() - 1) {
                        density << ',';
                    }
                }
                density << std::endl;
            }
            density.close();
        }
        flush(std::cout);
        flush(population);
        flush(moment1);
        flush(density);
    }
    population.close();
    moment1.close();

    density.open("density.csv");
    for (auto row : grid.cells) {
        for (int i = 0; i < (int)row.size(); ++i) {
            density << row[i].population;
            if (i < (int)row.size() - 1) {
                density << ',';
            }
        }
        density << std::endl;
    }
    density.close();

    return 0;
}
