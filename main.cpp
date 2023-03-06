#include "PolynomEvaluation.h"
#include <iomanip>

int main() {
    std::vector<double> roots = {1, 1, 1, 1, 1};

    Polynom<double, 5> p(polynomial_coeffs(roots));

    double x = 1;

    std::cout << "||" << calc_error(p, x) << "||"<< std::endl;

    bool result = fabs(horner(p, x)) < calc_error(p, x);

    std::cout << std::fixed << std::setprecision(16) << result;

}
