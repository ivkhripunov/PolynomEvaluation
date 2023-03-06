#include "../src/PolynomEvaluation.h"
#include <gtest/gtest.h>
#include <fstream>

TEST(POLYNOM_EVAL, TEST_1) {
    std::ofstream file;
    file.open("/home/ivankhripunov/CLionProjects/PolynomEvaluation/tests/log_1.txt");

    std::vector<double> roots = {1, 1, 1, 1, 1};
    Polynom<double, 5> p(polynomial_coeffs(roots));

    double step = 1e-6;
    boost::multiprecision::cpp_dec_float_100 boost_step(1e-6);

    for (std::size_t i = 0; i < 2e5; ++i) {
        boost::multiprecision::cpp_dec_float_100 boost_result("0");
        double comp_horner_result = 0, horner_result = 0, basic_result = 0;

        double x = 0.9000000 + step * i;
        boost::multiprecision::cpp_dec_float_100 boost_x = 0.9000000 + i * boost_step;

        for (std::size_t j = 0; j < 10; ++j) {

            basic_result += basic_evaluation(p, x);
            horner_result += horner(p, x);
            comp_horner_result += compensated_horner(p, x);
            boost_result += basic_evaluation(p, boost_x);

        }

        file << std::fixed << std::setprecision(16) << x << " " << (basic_result - boost_result) / 100 << " "
                  << (horner_result - boost_result) / 100 << " " << (comp_horner_result - boost_result) / 100 << std::endl;
    }

    file.close();

}

