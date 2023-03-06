#include "../src/PolynomEvaluation.h"
#include <gtest/gtest.h>


TEST(POLYNOM_EVAL, TEST_1) {

    boost::multiprecision::cpp_dec_float_100 boost_x("0.9993746556");
    boost::multiprecision::cpp_dec_float_100 boost_result("0");


    std::vector<double> roots = {1, 1, 1, 1, 1};
    double comp_horner_result, horner_result, basic_result;

    Polynom<double, 5> p(polynomial_coeffs(roots));

    double x = boost::numeric::converter<double, boost::multiprecision::cpp_dec_float_100>::convert(boost_x);

    std::cout << "||" << calc_error(p, x) << "||" << std::endl;


    basic_result = basic_evaluation(p, x);
    horner_result = horner(p, x);
    comp_horner_result = compensated_horner(p, x);
    boost_result = basic_evaluation(p, boost_x);

    std::cout << std::fixed << std::setprecision(16) << (basic_result - boost_result) << " "
              << horner_result - boost_result << " " << comp_horner_result - boost_result << std::endl;

}

