#include "PolynomEvaluation.h"
#include <iomanip>

int main() {
    SparsePolynom<double> p = {{0,  1},
                               };

    double result = compensated_horner(p, 1.0001);

    std::cout << std::fixed << std::setprecision(16) << result;
}
