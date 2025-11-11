cmake_minimum_required(VERSION 3.31)
project(FindingGenerators)
set(CMAKE_CXX_STANDARD 23)

set(Python3_EXECUTABLE "/Users/oliver/Documents/programming/FindingGenerators/.venv/bin/python3")
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

execute_process(
        COMMAND ${Python3_EXECUTABLE} -m pybind11 --cmakedir
        OUTPUT_VARIABLE PYBIND11_CMAKE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
)

set(pybind11_DIR ${PYBIND11_CMAKE_DIR})

message(STATUS "Using pybind11 found at: ${pybind11_DIR}")

add_subdirectory(pybind11)
pybind11_add_module(fg SHARED main.cpp)

set(INIT_PY "${CMAKE_SOURCE_DIR}/package/__init__.py")
if(NOT EXISTS ${INIT_PY})
    file(WRITE ${INIT_PY} "")
endif()

add_custom_command(
        TARGET fg
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_SOURCE_DIR}/package
        COMMAND ${CMAKE_COMMAND} -E copy
        $<TARGET_FILE:fg>
        ${CMAKE_SOURCE_DIR}/package/
)