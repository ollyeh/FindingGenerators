#include <iostream>
#include <pybind11/pybind11.h>

float func(float x) {
    return x*x;
}

PYBIND11_MODULE(fg, m) {
    m.doc() = "Python bindings for functions in Python";
}
