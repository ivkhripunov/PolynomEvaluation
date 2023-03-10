#include "../src/PolynomEvaluation.h"
#include <gtest/gtest.h>
#include <fstream>

TEST(POLYNOM_EVAL, TEST_1) {
    std::ofstream file;
    file.open("/home/ivankhripunov/CLionProjects/PolynomEvaluation/tests/log_1.txt");

    std::vector<double> roots = {1, 2, 3, 4, 5, 6, 7};
    std::vector<boost::multiprecision::cpp_dec_float_100> boost_roots = {1, 2, 3, 4, 5, 6, 7};
    Polynom<double, 7> p(polynomial_coeffs(roots));

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
            boost_result += basic_evaluation(boost_p, boost_x);

        }

        file << std::fixed << std::setprecision(16)
             << x << " " << boost_result / n << " "
             << basic_result / n << " "
             << horner_result / n << " "
             << comp_horner_result / n << " "
             << calc_error(p, x, 1.11e-16) << std::endl;

    }

    file.close();

}

TEST(POLYNOM_EVAL, TEST_2) {
    std::ofstream file;
    file.open("/home/ivankhripunov/CLionProjects/PolynomEvaluation/tests/log_1.txt");

    std::vector<double> roots = {1, 1, 1, 1, 1};
    std::vector<boost::multiprecision::cpp_dec_float_100> boost_roots = {1, 1, 1, 1, 1};
    Polynom<double, 5> p(polynomial_coeffs(roots));

    Polynom<boost::multiprecision::cpp_dec_float_100, 7> boost_p(polynomial_coeffs(boost_roots));

    double step = 1e-4;
    boost::multiprecision::cpp_dec_float_100 boost_step(1e-4);

    for (std::size_t i = 0; i < 6e1; ++i) {
        boost::multiprecision::cpp_dec_float_100 boost_result("0");
        double comp_horner_result = 0, horner_result = 0, basic_result = 0;

        double x = 0.998 + step * i;
        boost::multiprecision::cpp_dec_float_100 start("0.998");
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
             << 0 << std::endl;
    }

    file.close();

}

TEST(POLYNOM_EVAL, TEST_3) {
    std::ofstream file;
    file.open("/home/ivankhripunov/CLionProjects/PolynomEvaluation/tests/log_2.txt");

    std::vector<double> roots = {1.5346, -23957, 145.567, 384.4657, -0.00245, 6.35887, -70.5};
    std::vector<boost::multiprecision::cpp_dec_float_100> boost_roots = {1.5346, -23957, 145.567, 384.4657, -0.00245, 6.35887, -70.5};
    Polynom<double, 7> p(polynomial_coeffs(roots));

    Polynom<boost::multiprecision::cpp_dec_float_100, 7> boost_p(polynomial_coeffs(boost_roots));

    double step = 1e-4;
    boost::multiprecision::cpp_dec_float_100 boost_step(1e-4);

    for (std::size_t i = 0; i < 1e2; ++i) {
        boost::multiprecision::cpp_dec_float_100 boost_result("0");
        double comp_horner_result = 0, horner_result = 0, basic_result = 0;

        double x = 1.53 + step * i;
        boost::multiprecision::cpp_dec_float_100 start("1.53");
        boost::multiprecision::cpp_dec_float_100 boost_x = start + i * boost_step;



        basic_result = basic_evaluation(p, x);
        horner_result = horner(p, x);
        comp_horner_result = compensated_horner(p, x);
        boost_result = horner(boost_p, boost_x);


        file << std::fixed << std::setprecision(16)
             << x << " " << boost_result << " "
             << horner_result << " "
             << comp_horner_result << std::endl;
    }

    file.close();

}

const std::size_t deg = 40;
const double x = 1.333;

TEST(POLYNOM_EVAL, TEST_4) {
    std::ofstream out;
    std::ifstream in;
    std::string line;
    double tmp;
    out.open("/home/ivankhripunov/CLionProjects/PolynomEvaluation/tests/log_1.txt", std::ios::app);
    in.open("/home/ivankhripunov/CLionProjects/PolynomEvaluation/tests/poly_coeff.txt");

    std::vector<double> polynom_coeff;

    for (std::size_t i = 0; i < deg + 1; ++i) {
        in >> tmp;
        polynom_coeff.push_back(tmp);
    }

    Polynom<double, deg> p(polynom_coeff);

    double comp_horner_result, horner_result;

    horner_result = horner(p, x);
    comp_horner_result = compensated_horner(p, x);

    out << std::fixed << std::setprecision(23) << horner_result << " " << comp_horner_result << " " << deg << std::endl;


    out.close();
    in.close();

}
