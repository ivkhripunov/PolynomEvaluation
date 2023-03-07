#ifndef POLYNOMEVALUATION_POLYNOMEVALUATION_H
#define POLYNOMEVALUATION_POLYNOMEVALUATION_H

#include "Polynom.h"
#include <boost/multiprecision/cpp_dec_float.hpp>

template<typename Type>
struct ReturnStruct {
    Type result;
    Type error;
};

/**
 * Error-free transformation of the sum of 2 floating point numbers
 * @tparam Type floating point type
 * @param a floating point number
 * @param b floating point number
 * @return struct: sum of numbers and error
 */
template<typename Type>
ReturnStruct<Type> two_sum(const Type &a,
                           const Type &b) {

    ReturnStruct<Type> out;

    out.result = a + b;
    Type tmp = out.result - a;
    out.error = a - (out.result - tmp) + (b - tmp);

    return out;
}

/**
 * Error-free transformation of the product of to floating point numbers with Fused Multiply and add (FMA)
 * @tparam Type floating point type
 * @param a floating point number
 * @param b floating point number
 * @return struct: result of a * b and error
 */
template<typename Type>
ReturnStruct<Type> two_product_fma(const Type &a,
                                   const Type &b) {
    ReturnStruct<Type> out;
    out.result = a * b;
    out.error = std::fma(a, b, -out.result);

    return out;
}

/**
 * Horner scheme
 * @tparam Type floating point type
 * @tparam N polynom degree
 * @param polynom polynom with FP coeffs
 * @param x value for polynom calculation
 * @return polynom value in point x
 */
template<typename Type, std::size_t N>
Type horner(const Polynom<Type, N> &polynom, const Type &x) {

    Type sum = polynom[polynom.get_degree()], p;

    for (long long i = polynom.get_degree() - 1; i >= 0; i--) {
        p = sum * x;
        sum = p + polynom[i];
    }

    return sum;
}

/**
 * Compensated Horner Scheme
 * @tparam Type floating point type
 * @tparam N polynom degree
 * @param polynom polynom with FP coeffs
 * @param x value for polynom calculation
 * @return polynom value in point x
 */
template<typename Type, std::size_t N>
Type
compensated_horner(const Polynom<Type, N> &polynom, const Type &x) {

    Polynom<Type, N - 1> polynom_pi, polynom_sigma;

    ReturnStruct<Type> p, s;
    s.result = polynom[polynom.get_degree()];

    for (long long i = polynom.get_degree() - 1; i >= 0; i--) {

        p = two_product_fma(s.result, x);
        s = two_sum(p.result, polynom[i]);

        polynom_pi[i] = p.error;
        polynom_sigma[i] = s.error;
    }

    return horner(polynom, x) + horner(polynom_pi + polynom_sigma, x);
}

/**
 * Calculating polynom coeffs from given roots
 * @tparam Type floating point type
 * @param roots polynom roots
 * @return vector of coeffs, vector[0] is zero degree
 */
template<typename Type>
std::vector<Type> polynomial_coeffs(const std::vector<Type> &roots) {

    std::vector<Type> coeffs(roots.size() + 1, 0);
    coeffs[0] = static_cast<Type>(1.);

    for (std::size_t i = 0; i < roots.size(); ++i) {
        for (std::size_t j = roots.size(); j >= 1; j--) {
            coeffs[j] = coeffs[j - 1] - roots[i] * coeffs[j];
        }

        coeffs[0] *= -roots[i];
    }

    return coeffs;

}

/**
 * Calculating condition number
 * @tparam Type floating point type
 * @tparam N polynom degree
 * @param polynom polynom with FP coeffs
 * @param x value for polynom calculation
 * @return condition number
 */
template<typename Type, std::size_t N>
Type calc_condition_number(const Polynom<Type, N> &polynom, const Type &x) {
    Type sum_1 = abs(horner(polynom, x));
    Polynom<Type, N> abs_polynom(polynom);

    for (std::size_t i = 0; i < N + 1; ++i) abs_polynom[i] = fabs(polynom[i]);

    Type sum_2 = horner(abs_polynom, fabs(x));

    return sum_2 / sum_1;
}

/**
 * Auxilary function for error calculation
 */
template<typename Type>
Type calc_gamma(const std::size_t &n, const Type &u) {
    return n * u / (1 - n * u);
}

/**
 * Auxiliary function
 */
template<typename Type>
Type power(const Type &x, const std::size_t &degree) {
    Type result = 1;
    for (std::size_t i = 0; i < degree; ++i) result *= x;

    return result;
}

/**
 * Basic way of evaluating polynom
 * @tparam Type floating point type
 * @tparam N polynom degree
 * @param polynom polynom with FP coeffs
 * @param x value for polynom calculation
 * @return polynom value in point x
 */
template<typename Type, std::size_t N>
Type basic_evaluation(const Polynom<Type, N> &polynom, const Type &x) {
    Type sum = static_cast<Type>(0);

    for (std::size_t i = 0; i < polynom.get_degree() + 1; ++i) sum += polynom[i] * power(x, i);

    return sum;
}

#endif //POLYNOMEVALUATION_POLYNOMEVALUATION_H
