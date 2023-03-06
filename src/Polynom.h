#ifndef POLYNOMEVALUATION_POLYNOM_H
#define POLYNOMEVALUATION_POLYNOM_H

#include <array>
#include <vector>
#include <iostream>
#include <initializer_list>
#include <cmath>
#include <iterator>

template<typename Type, std::size_t N>
class Polynom {

protected:
    std::array<Type, N + 1> data_;
public:
    constexpr Polynom() {
        data_.fill(static_cast<Type>(0));
    }

    constexpr Polynom(const Polynom<Type, N> &other) {
        for (auto i = 0; i < N + 1; ++i) data_[i] = *(other.data_.begin() + i);
    }

    constexpr Polynom(const std::initializer_list<Type> &initializer) {
        for (auto i = 0; i < initializer.size(); ++i) data_[i] = *(initializer.begin() + i);
    }

    constexpr Polynom(const std::vector<Type> &coeffs) {
        for (std::size_t i = 0; i < data_.size(); ++i) data_[i] = coeffs[i];
    }

    const Type &operator[](std::size_t i) const {
        return data_[i];
    }

    Type &operator[](std::size_t i) {
        return data_[i];
    }

    [[nodiscard]] std::size_t get_degree() const {
        return N;
    }

    Polynom<Type, N> operator+=(const Polynom<Type, N> &other) {
        for (auto i = 0; i < N; ++i) data_[i] += other[i];
        return *this;
    }

    Polynom<Type, N> operator+(const Polynom<Type, N> &other) const {
        Polynom<Type, N> result(*this);

        result += other;

        return result;
    }

    friend std::ostream &operator<<(std::ostream &out, const Polynom<Type, N> &vec) {
        for (const auto &element: vec.data_) out << element << " ";
        out << std::endl;

        return out;
    }
};

#endif //POLYNOMEVALUATION_POLYNOM_H
