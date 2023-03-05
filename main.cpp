#include "PolynomEvaluation.h"

int main() {
    SparsePolynom<double> p = {{0, 2},
                               {1, -1},
                               {2, 5}};

    double result = horner_scheme(p, 3.);

    std::cout << result;
}
