#include "../src/PolynomEvaluation.h"
#include <gtest/gtest.h>
#include <fstream>

TEST(POLYNOM_EVAL, TEST_1) {
    std::ofstream file;
    file.open("/home/ivankhripunov/CLionProjects/PolynomEvaluation/tests/log_1.txt");

    std::vector<double> roots = {1, 2, 3, 4, 5, 6, 7};
    std::vector<boost::multiprecision::cpp_dec_float_100> boost_roots = {1, 2, 3, 4, 5, 6, 7};
    Polynom<double, 7> p(polynomial_coeffs(roots));
    std::cout << p;
    Polynom<boost::multiprecision::cpp_dec_float_100, 7> boost_p(polynomial_coeffs(boost_roots));

    double step = 1e-4;
    boost::multiprecision::cpp_dec_float_100 boost_step(1e-4);

    for (std::size_t i = 0; i < 2e2; ++i) {
        boost::multiprecision::cpp_dec_float_100 boost_result("0");
        double comp_horner_result = 0, horner_result = 0, basic_result = 0;

        double x = 1.99 + step * i;
        boost::multiprecision::cpp_dec_float_100 start("1.99");
        boost::multiprecision::cpp_dec_float_100 boost_x = start + i * boost_step;

        std::size_t n = 100;

        for (std::size_t j = 0; j < n; ++j) {

            basic_result += basic_evaluation(p, x);
            horner_result += horner(p, x);
            comp_horner_result += compensated_horner(p, x);
            boost_result += horner(boost_p, boost_x);

        }

        file << std::fixed << std::setprecision(16)
             << x << " " << boost_result / n << " "
             << basic_result / n << " "
             << horner_result / n << " "
             << comp_horner_result / n << " "
             << calc_error(p, x) << std::endl;

        //ASSERT_TRUE(abs(comp_horner_result - boost_result) / n < calc_error(p, x));
    }

    file.close();

}

