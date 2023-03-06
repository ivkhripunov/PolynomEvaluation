#ifndef POLYNOMEVALUATION_POLYNOMEVALUATION_H
#define POLYNOMEVALUATION_POLYNOMEVALUATION_H

#include "Polynom.h"

template<typename Type>
Type two_sum(Type &error, const Type &a,
             const Type &b) {
    /// Функция для вычисления суммы a + b и ошибки на шаге

    Type sum = a + b; //вычисляем сумму напрямую
    Type b_virtual = sum - a; // часть "b" в "sum" : =b if a<b
    Type a_virtual = sum - b_virtual; // часть "a" в "sum" : =a if a>b

    Type b_round_error = b - b_virtual; //ошибка округления "b" : =0 if a<b
    Type a_round_error = a - a_virtual; //ошибка округления "a" : =0 if a>b

    error += a_round_error + b_round_error; //к итоговой ошибке прибавляем ошибку на шаге

    return sum; //возвращаем сумму
}

template<typename Type>
Type split(Type &error, const Type &a) {
    int s = 53 / 2;
    Type factor = pow(2, s) + 1;
    Type c = factor * a;
    Type result = c - (c - a);
    error = a - result;

    return result;
}

template<typename Type>
Type two_product_fma(Type &error, const Type &a,
                     const Type &b) {
    Type result = a * b;
    error = std::fma(a, b, -result);
}

template<typename Type, std::size_t N>
Type horner(const Polynom<Type, N> &polynom, const Type &x) {

    Type sum = polynom[polynom.get_degree()];

    for (long long i = polynom.get_degree() - 1; i >= 0; i--) {
        sum = polynom[i] + sum * x;
    }

    return sum;
}

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

template<typename Type, std::size_t N>
Type compensated_horner(const Polynom<Type, N> &polynom, const Type &x) {

    Polynom<Type, N - 1> p_pi, p_sigma;

    Type h = EFT_horner(polynom, x, p_pi, p_sigma);

    Polynom<Type, N - 1> poly = p_pi + p_sigma;

    std::cout << poly;

    Type c = horner(poly, x);
    Type res = h + c;

    return res;
}

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

template<typename Type, std::size_t N>
Type calc_condition_number(const Polynom<Type, N> &polynom, const Type &x) {
    Type sum_1 = abs(horner(polynom, x));
    Polynom<Type, N> abs_polynom(polynom);

    for (std::size_t i = 0; i < N + 1; ++i) abs_polynom[i] = fabs(polynom[i]);

    Type sum_2 = horner(abs_polynom, fabs(x));

    return sum_2 / sum_1;
}

template<typename Type>
Type calc_gamma(const std::size_t &n, const Type &u) {
    return n * u / (1 - n * u);
}

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

#endif //POLYNOMEVALUATION_POLYNOMEVALUATION_H
