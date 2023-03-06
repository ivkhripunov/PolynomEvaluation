#include "../src/PolynomEvaluation.h"
#include <gtest/gtest.h>


TEST(POLYNOM_EVAL, TEST_1) {
    std::vector<double> roots = {1, 1, 1, 1, 1};

    Polynom<double, 5> p(polynomial_coeffs(roots));

    double x = 0.999999999;

    std::cout << "||" << calc_error(p, x) << "||"<< std::endl;

    std::cout << basic_evaluation(p, x) << " "<< horner(p, x) <<" " << compensated_horner(p, x);

    ASSERT_NEAR(basic_evaluation(p, x), -3.3306690738754696e-16, calc_error(p, x));

}

