#include "PolynomEvaluation.h"
#include <iomanip>

int main() {
    SparsePolynom<double> p = {{0,  1},
                               {1,  -1},
                               {2,  1},
                               {3,  -1},
                               {4,  1},
                               {5,  -1},
                               {6,  1},
                               {7,  -1},
                               {8,  1},
                               {9,  -1},
                               {10, 1}};

    double result = compensated_horner(p, 1.0001);

    std::cout << std::fixed << std::setprecision(16) << result;
}
