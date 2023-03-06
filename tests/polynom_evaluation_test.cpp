#include "../src/PolynomEvaluation.h"
#include <gtest/gtest.h>


TEST(POLYNOM_EVAL, TEST_1) {
    std::vector<double> roots = {1, 1, 1, 1, 1};

    Polynom<double, 5> p(polynomial_coeffs(roots));

    double x = 2;

    std::cout << "||" << calc_error(p, x) << "||"<< std::endl;

    ASSERT_NEAR(compensated_horner(p, x), 1., calc_error(p, x));

}

