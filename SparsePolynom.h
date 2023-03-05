#ifndef POLYNOMEVALUATION_SPARSEPOLYNOM_H
#define POLYNOMEVALUATION_SPARSEPOLYNOM_H

#include <map>
#include <iostream>
#include <initializer_list>
#include <algorithm>
#include <cmath>

template<typename Type>
class SparsePolynom {

private:

    std::map<std::size_t, Type> data_; //пары степень - коэффициент
    Type zero_ = static_cast<Type>(0);

public:

    SparsePolynom() = default;

    SparsePolynom(const SparsePolynom<Type> &other) {
        for (const auto& pair : data_) data_.insert(pair);
    }

    SparsePolynom(const std::initializer_list<std::pair<std::size_t, Type>> &initializer) {
        for (const auto &pair: initializer) data_.insert(pair);
    }

    const Type &operator[](const std::size_t &i) const {
        if (data_.contains(i)) return data_.at(i);
        else return zero_;
    }

    bool contains_degree(const std::size_t &degree) const {
        return data_.contains(degree);
    }


    [[nodiscard]] std::size_t get_degree() const {
        if (data_.empty()) return 0;
        return (*std::find_if(std::rbegin(data_), std::rend(data_), [](auto &element) {
            return std::fabs(element.second - 0) > std::numeric_limits<Type>::epsilon();
        })).first;
    }

    void add(const std::pair<std::size_t, Type> &pair) {
        data_.insert(pair);
    }

    SparsePolynom<Type> operator+=(const SparsePolynom<Type> &other) {
        for (const auto &pair: other.data_) {
            if (data_.contains(pair.first)) data_.at(pair.first) += pair.second;
            else data_.insert(pair);
        }
        return *this;
    }

    SparsePolynom<Type> operator+(const SparsePolynom<Type> &other) const {
        SparsePolynom<Type> result(*this);
        result += other;

        return result;
    }


    friend std::ostream &operator<<(std::ostream &out, const SparsePolynom<Type> &polynom) {
        for (const auto &pair: polynom.data_) std::cout << "[" << pair.first << ", " << pair.second << "]";
        std::cout << std::endl;

        return out;
    }

};

#endif //POLYNOMEVALUATION_SPARSEPOLYNOM_H
