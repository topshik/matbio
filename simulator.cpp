#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

const double birth_rate = 0.4;
const double death_rate = 0.2;

class Cell {
public:
    int population;
    double x, y;
    int row, column;
    Cell() {};
    Cell(int Population, double X, double Y)
    : population(Population), x(X), y(Y) {};
};

class Grid {
public:
    int width, height;
    int n;
    int total_population;
    std::vector<std::vector<Cell> > cells;
    Grid() {};
    Grid(int Width, int Height, int N)
    : width(Width), height(Height), n(N) {
        std::vector<Cell> init_cells;
        Cell init_cell;
        total_population = (n + 1) * (n + 1);
        for (size_t i = 0; i != n; ++i) {
            cells.push_back(init_cells);
            for (size_t j = 0; j != n; ++j) {          // grid is inited with 1 unit in every cell
                init_cell.population = 1;
                init_cell.x = (width / n) * i;
                init_cell.y = (height / n) * j;
                init_cell.row = i;
                init_cell.column = j;
                cells[i].push_back(init_cell);
            }
        }
    };

    Cell cell(int i, int j) {
        return cells[i][j];
    }

    int cell_population(int i, int j) {
        return cells[i][j].population;
    }
};

double distance(const Cell &from, const Cell &to) {
    return std::sqrt(std::pow((from.x - to.x), 2) + std::pow((from.y - to.y), 2));
}

int max_distance = 2.5;
std::vector<Cell> neighbour_birth_influence(const Grid &grid, const Cell &cell) {
    std::vector<Cell> result;
    int border_x = max_distance * grid.n / grid.width;  // Max distance from cell by x
    int border_y = max_distance * grid.n / grid.height;  // Max distance from cell by y
    for (size_t i = (cell.row - border_x < 0 ? 0 : cell.row - border_x); i <
                    (cell.row + border_x > grid.width ? 0 : cell.row + border_x); ++i) {
        for (size_t j = (cell.column - border_y < 0 ? 0 : cell.column - border_y); i <
                    (cell.column + border_y > grid.height ? 0 : cell.column + border_y); ++i) {
            if (distance(cell, grid.cells[i][j]) < max_distance) {
                result.push_back(grid.cells[i][j]);
            }
        }
    }
    return result;
};

void iteration(Grid &grid) {
    std::vector<std::vector<double> > nobirth_matrix(grid.height, std::vector<double>(grid.width, 1));

    /* Counting (no)birth probabilities */
    for (auto row : grid.cells) {
        for (auto cell : row) {
            std::vector<Cell> neighbours = neighbour_birth_influence(grid, cell);
            for (auto neighbour : neighbours) {
                double birth_prob = birth_rate * std::exp(-std::pow(distance(neighbour, cell), 2));
                nobirth_matrix[neighbour.row][neighbour.column] *= (1 - birth_prob);
            }
        }
    }

    /* Giving birth and killing*/
    for (size_t i = 0; i < nobirth_matrix.size(); ++i) {
        for (size_t j = 0; j < nobirth_matrix[0].size(); ++j) {
            if ((double)rand() / RAND_MAX > nobirth_matrix[i][j]) {
                ++grid.cells[i][j].population;
                ++grid.total_population;
            }

            /* Use binominal_distribution instead! */
            int died = 0;
            for (size_t k = 0; k < grid.cells[i][j].population; ++k) {
                if ((double)rand() / RAND_MAX < death_rate) {
                    ++died;
                }
            }
            grid.cells[i][j].population -= died;
            grid.total_population -= died;
        }
    }
}

int main() {
    Grid grid(10, 10, 10000);

    size_t iterations = 10000;
    for (size_t i = 0; i < iterations; ++i) {
        iteration(grid);
        std::cout << i << "\t" << grid.total_population << std::endl;
    }
}
