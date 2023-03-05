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

public:

    SparsePolynom() = default;

    SparsePolynom(const SparsePolynom<Type> &other) {
        data_ = other.data_;
    }

    SparsePolynom(const std::initializer_list<std::pair<std::size_t, Type>> &initializer) {
        for (const auto &pair: initializer) data_.insert(pair);
    }

    const Type &operator[](const std::size_t &i) const {
        if (data_.contains(i)) return data_.at(i);
        else return static_cast<Type>(0);
    }

    bool contains_degree(const std::size_t &degree) const {
        return data_.contains(degree);
    }


    [[nodiscard]] std::size_t get_degree() const {
        return (*std::find_if(std::rbegin(data_), std::rend(data_), [](auto &element) {
            return std::fabs(element.second - 0) > std::numeric_limits<double>::epsilon();
        })).first;
    }


    friend std::ostream &operator<<(std::ostream &out, const SparsePolynom<Type> &polynom) {
        for (const auto &pair: polynom.data_) std::cout << "[" << pair.first << ", " << pair.second << "]";
        std::cout << std::endl;

        return out;
    }

};

#endif //POLYNOMEVALUATION_SPARSEPOLYNOM_H
