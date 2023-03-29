#ifndef POLYNOMEVALUATION_POLYNOM_H
#define POLYNOMEVALUATION_POLYNOM_H

#include <array>
#include <cmath>

using indexType = std::size_t;
using scalar = double;

/*** Функцию для подсчета ***/

namespace Containers {
    template<typename T, indexType N>
    using array = std::array<T, N>;
}

template<typename T, indexType N>
class Polynom {

private:
    Containers::array<T, N + 1> data_;

public:

    constexpr Polynom() {
        data_.fill(static_cast<T>(0));
    }

    constexpr Polynom(const Containers::array<T, N + 1> &coeffs) {
        data_ = coeffs;
    }

    const T &operator[](const indexType &i) const {
        return data_[i];
    }

    T &operator[](const indexType &i) {
        return data_[i];
    }


    Polynom<T, N> operator+(const Polynom<T, N> &other) const {
        Polynom<T, N> result(*this);

        for (auto i = 0; i < N + 1; ++i) {
            result[i] += other[i];
        }

        return result;
    }
};

template<typename T>
struct ReturnStruct {
    T result;
    T error;
};

/**
 * Error-free transformation of the sum of 2 floating point numbers
 * @tparam T floating point type
 * @param a floating point number
 * @param b floating point number
 * @return struct: sum of numbers and error
 */
template<typename T>
ReturnStruct<T> TwoSum(const T &a, const T &b) {

    ReturnStruct<T> out;

    out.result = a + b;
    const T b_virtual = out.result - a;
    const T a_virtual = out.result - b_virtual;
    const T b_roundoff = b - b_virtual;
    const T a_roundoff = a - a_virtual;
    out.error = a_roundoff + b_roundoff;

    return out;
}

/**
 * Error-free transformation of the product of to floating point numbers with Fused Multiply and add (FMA)
 * @tparam T floating point type
 * @param a floating point number
 * @param b floating point number
 * @return struct: result of a * b and error
 */
template<typename T>
ReturnStruct<T> TwoProductFMA(const T &a,
                              const T &b) {
    ReturnStruct<T> out;

    out.result = a * b;
    out.error = std::fma(a, b, -out.result);

    return out;
}

/**
 * Horner scheme
 * @tparam T floating point type
 * @tparam N polynom degree
 * @param polynom polynom with FP coeffs
 * @param x value for polynom calculation
 * @return polynom value in point x
 */
template<typename T, indexType N>
T Horner(const Polynom<T, N> &polynom, const T &x) {

    T sum = polynom[N];

    for (indexType i = N; i >= 1; i--) {
        sum = sum * x + polynom[i - 1];
    }

    return sum;
}

/**
 * Compensated Horner Scheme
 * @tparam T floating point type
 * @tparam N polynom degree
 * @param polynom polynom with FP coeffs
 * @param x value for polynom calculation
 * @return polynom value in point x
 */
template<typename T, indexType N>
T CompensatedHorner(const Polynom<T, N> &polynom, const T &x) {

    Polynom<T, N - 1> polynom_pi, polynom_sigma;

    ReturnStruct<T> p, s;
    s.result = polynom[N];

    for (indexType i = N; i >= 1; i--) {

        p = TwoProductFMA(s.result, x);
        s = TwoSum(p.result, polynom[i - 1]);

        polynom_pi[i - 1] = p.error;
        polynom_sigma[i - 1] = s.error;
    }

    return s.result + Horner(polynom_pi + polynom_sigma, x);
}

#endif //POLYNOMEVALUATION_POLYNOM_H
