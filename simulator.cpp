#include <iostream>
#include <vector>

const double birth_rate = 0.4;
const double death_rate = 0.2;

class Cell {
public:
    int population;
    double x, y;
    std::vector<double> death_rates;
    Cell() {};
    Cell(int Population, double X, double Y)
    : population(Population), x(X), y(Y) {
    death_rates = std::vector<double> (population, death_rate);
    };
};

class Grid {
public:
    int length_x, length_y;
    int n;
    int total_population;
    std::vector<std::vector<Cell> > cells;
    Grid() {};
    Grid(int Length_x, int Length_y, int N)
    : length_x(Length_x), length_y(Length_y), n(N) {
        std::vector<Cell> init_cells;
        Cell init_cell;
        for (size_t i = 0; i != n; ++i) {
            cells.push_back(init_cells);
            for (size_t j = 0; j != n; ++j) {
                init_cell.population = 1;
                init_cell.x = (length_x / n) * i;
                init_cell.y = (length_y / n) * j;
                cells[i].push_back(init_cell);
            }
        }
    };

    Cell cell(int i, int j)
    {
        return cells[i][j];
    }

    int cell_population(int i, int j)
    {
        return cells[i][j].population;
    }
};
/*
void iteration(Grid grid) {
    for ()
}
*/
int main() {

}
