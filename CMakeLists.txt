cmake_minimum_required(VERSION 3.31)
project(FindingGenerators)

set(CMAKE_CXX_STANDARD 23)

add_subdirectory(pybind11)
pybind11_add_module(FindingGenerators main.cpp)