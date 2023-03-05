#ifndef POLYNOMEVALUATION_POLYNOMEVALUATION_H
#define POLYNOMEVALUATION_POLYNOMEVALUATION_H

#include "SparsePolynom.h"

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
Type two_product(const Type &a,
                 const Type &b) {
    return a * b;
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
Type two_product(Type &error, const Type &a,
                 const Type &b) {
    Type result = a * b;

    Type a_, error_a_;
    Type b_, error_b_;

    a_ = split(error_a_, a);
    b_ = split(error_b_, b);

    error = error_a_ * error_b_ - (((result - a_ * b_) - error_a_ * b_) - a_ * error_b_);

    return result;
}

template<typename Type>
Type two_product_fma(Type &error, const Type &a,
                     const Type &b) {
    Type result = a * b;
    error = std::fma(a, b, -result);
}

template<typename Type>
Type horner(const SparsePolynom<Type> &polynom, const Type &x) {

    Type sum = polynom[polynom.get_degree()];

    for (long long i = polynom.get_degree() - 1; i >= 0; i--) {
        sum = polynom[i] + sum * x;
    }

    return sum;
}

template<typename Type>
Type
EFT_horner(const SparsePolynom<Type> &polynom, const Type &x, SparsePolynom<Type> &p_pi, SparsePolynom<Type> &p_sigma) {

    Type p_i, pi_i, s_i = polynom[polynom.get_degree()], sigma_i;

    for (long long i = polynom.get_degree() - 1; i >= 0; i--) {

        p_i = two_product_fma(pi_i, s_i, x);
        s_i = two_sum(sigma_i, p_i, polynom[i]);

        p_pi.add({i, pi_i});
        p_sigma.add({i, sigma_i});
    }

    return horner(polynom, x);
}

template<typename Type>
Type compensated_horner(const SparsePolynom<Type> &polynom, const Type &x) {

    SparsePolynom<Type> p_pi, p_sigma;

    Type h = EFT_horner(polynom, x, p_pi, p_sigma);

    SparsePolynom<Type> poly = p_pi + p_sigma;

    Type c = horner(poly, x);
    Type res = h + c;

    return res;
}


#endif //POLYNOMEVALUATION_POLYNOMEVALUATION_H
