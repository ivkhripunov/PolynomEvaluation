#ifndef POLYNOMEVALUATION_POLYNOMEVALUATION_H
#define POLYNOMEVALUATION_POLYNOMEVALUATION_H

#include "Polynom.h"

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
ReturnStruct<Type> two_sum(const Type &a, const Type &b) {

    ReturnStruct<Type> out;

    out.result = a + b;
    const Type b_virtual = out.result - a;
    const Type a_virtual = out.result - b_virtual;
    const Type b_roundoff = b - b_virtual;
    const Type a_roundoff = a - a_virtual;
    out.error = a_roundoff + b_roundoff;

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

    Type sum = polynom[polynom.get_degree()];

    for (long long i = polynom.get_degree() - 1; i >= 0; i--) {
        sum = sum * x + polynom[i];
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
Type compensated_horner(const Polynom<Type, N> &polynom, const Type &x) {

    Polynom<Type, N - 1> polynom_pi, polynom_sigma;

    ReturnStruct<Type> p, s;
    s.result = polynom[N];

    for (long long i = N - 1; i >= 0; i--) {

        p = two_product_fma(s.result, x);
        s = two_sum(p.result, polynom[i]);

        polynom_pi[i] = p.error;
        polynom_sigma[i] = s.error;
    }

    return s.result + horner(polynom_pi + polynom_sigma, x);
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

    return compensated_horner(polynom.get_abs(), fabs(x)) / fabs(compensated_horner(polynom, x));
}

/**
 * Auxilary function for error calculation
 */
template<typename Type>
Type calc_gamma(const std::size_t &n, const Type &unit_roundoff) {
    return n * unit_roundoff / (1 - n * unit_roundoff);
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

template<typename Type, std::size_t N>
Type calc_error(const Polynom<Type, N> &polynom, const Type &x, const Type &unit_roundoff) {

    Polynom<Type, N - 1> polynom_pi, polynom_sigma;

    ReturnStruct<Type> p, s;
    s.result = polynom[polynom.get_degree()];

    for (long long i = polynom.get_degree() - 1; i >= 0; i--) {

        p = two_product_fma(s.result, x);
        s = two_sum(p.result, polynom[i]);

        polynom_pi[i] = p.error;
        polynom_sigma[i] = s.error;
    }

    return unit_roundoff * fabs(compensated_horner(polynom, x)) +
           gamma(4 * polynom.get_degree() + 2) *
           horner(polynom_pi.get_abs() + polynom_sigma.get_abs(), fabs(x)) +
           2 * unit_roundoff * unit_roundoff *
           fabs(compensated_horner(polynom, x));
}


#endif //POLYNOMEVALUATION_POLYNOMEVALUATION_H
