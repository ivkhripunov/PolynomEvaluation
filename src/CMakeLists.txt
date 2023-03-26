cmake_minimum_required(VERSION 3.24)
project(PolynomEvaluation)

set(CMAKE_CXX_STANDARD 20)

file(GLOB_RECURSE source *.h *.cpp *.hpp)

add_library(PolynomEvaluation INTERFACE ${source} PolynomEvaluation.h)