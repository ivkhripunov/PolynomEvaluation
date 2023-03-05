#ifndef POLYNOMEVALUATION_POLYNOMEVALUATION_H
#define POLYNOMEVALUATION_POLYNOMEVALUATION_H

#include "SparsePolynom.h"

template <typename Type>
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

template <typename Type>
Type two_product(const Type &a,
                    const Type &b) {
    return a * b;
}

template<typename Type>
Type horner_scheme(const SparsePolynom<Type> &polynom, const Type &x) {

    Type sum = polynom[polynom.get_degree()];

    for (long long i = polynom.get_degree() - 1; i >= 0; i--) {
        sum = polynom[i] + sum * x;
    }

    return sum;
}


#endif //POLYNOMEVALUATION_POLYNOMEVALUATION_H
