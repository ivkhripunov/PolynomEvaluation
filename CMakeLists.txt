cmake_minimum_required(VERSION 3.24)
project(PolynomEvaluation)

set(CMAKE_CXX_STANDARD 20)

add_executable(PolynomEvaluation main.cpp PolynomEvaluation.h SparsePolynom.h)
