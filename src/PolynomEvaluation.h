#ifndef POLYNOMEVALUATION_POLYNOMEVALUATION_H
#define POLYNOMEVALUATION_POLYNOMEVALUATION_H

#include "Polynom.h"
#include <boost/multiprecision/cpp_dec_float.hpp>

/**
 * Error-free transformation of the sum of 2 floating point numbers
 * @tparam Type floating point type
 * @param error operational error, is changed in func
 * @param a floating point number
 * @param b floating point number
 * @return sum of numbers
 */
template<typename Type>
Type two_sum(Type &error, const Type &a,
             const Type &b) {

    Type sum = a + b;
    Type tmp = sum - a;
    error = a - (sum - tmp) + (b - tmp);

    return sum; //возвращаем сумму
}

/**
 * Error-free transformation of the product of to floating point numbers with Fused Multiply and add (FMA)
 * @tparam Type floating point type
 * @param error operational error, is changed in func
 * @param a floating point number
 * @param b floating point number
 * @return result of a * b
 */
template<typename Type>
Type two_product_fma(Type &error, const Type &a,
                     const Type &b) {
    Type result = a * b;
    error = std::fma(a, b, -result);

    return result;
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
 * Error free transformation for the Horner scheme
 * @tparam Type floating point type
 * @tparam N polynom degree
 * @param polynom polynom with FP coeffs
 * @param x value for polynom calculation
 * @param p_pi error polynom
 * @param p_sigma error polynom
 * @return polynom evaluated via Horner scheme
 */
template<typename Type, std::size_t N>
Type
EFT_horner(const Polynom<Type, N> &polynom, const Type &x, Polynom<Type, N - 1> &p_pi, Polynom<Type, N - 1> &p_sigma) {

    Type p_i, pi_i, s_i = polynom[polynom.get_degree()], sigma_i;

    for (long long i = polynom.get_degree() - 1; i >= 0; i--) {

        p_i = two_product_fma(pi_i, s_i, x);
        s_i = two_sum(sigma_i, p_i, polynom[i]);

        p_pi[i] = pi_i;
        p_sigma[i] = sigma_i;
    }

    return horner(polynom, x);
}

/**
 * Compensated Horner scheme
 * @tparam Type floating point type
 * @tparam N polynom degree
 * @param polynom polynom with FP coeffs
 * @param x value for polynom calculation
 * @return polynom value in point x
 */
template<typename Type, std::size_t N>
Type compensated_horner(const Polynom<Type, N> &polynom, const Type &x) {

    Polynom<Type, N - 1> p_pi, p_sigma;

    Type h = EFT_horner(polynom, x, p_pi, p_sigma);

    Type c = horner(p_pi + p_sigma, x);

    return h + c;
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
 * Function for calculating maximum error of Horner scheme
 * u is chosen for double precision!!! Do not use for other types without setup!
 * @tparam Type floating point type
 * @tparam N polynom degree
 * @param polynom polynom with FP coeffs
 * @param x value for polynom calculation
 * @return max error
 */
template<typename Type, std::size_t N>
Type calc_error(const Polynom<Type, N> &polynom, const Type &x) {
    Type abs_res = fabs(horner(polynom, x));

    double u = 1.11e-16;

    Polynom<Type, N - 1> p_pi, p_sigma;

    EFT_horner(polynom, x, p_pi, p_sigma);

    Polynom<Type, N - 1> abs_p_pi(p_pi), abs_p_sigma(p_sigma);

    for (std::size_t i = 0; i < abs_p_pi.get_degree() + 1; ++i) {
        abs_p_pi[i] = fabs(p_pi[i]);
        abs_p_sigma[i] = fabs(p_sigma[i]);
    }

    Type error = u * abs_res + (gamma(4 * polynom.get_degree() + 2) * horner(abs_p_pi + abs_p_sigma, fabs(x)) +
                                2 * u * u * abs_res);

    return error;
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
