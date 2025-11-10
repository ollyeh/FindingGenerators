#include <iostream>
#include <pybind11/pybind11.h>

float func(float x) {
    return x*x;
}