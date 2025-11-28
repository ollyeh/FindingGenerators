cmake_minimum_required(VERSION 3.31)
project(FindingGenerators)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_OSX_DEPLOYMENT_TARGET 14.0)

# easy3d v2.6.1     https://github.com/LiangliangNan/Easy3D
# pybind11 v3.0.1   https://github.com/pybind/pybind11
# spdlog v1.16.0    https://github.com/gabime/spdlog
# nlohmann v3.12.0  https://github.com/nlohmann/json

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

# Easy3D
set(EASY3D_DIR "${CMAKE_SOURCE_DIR}/Easy3D")
message(STATUS "Using Easy3D found at: ${EASY3D_DIR}")
add_subdirectory(${EASY3D_DIR})

pybind11_add_module(finding_generators SHARED main.cpp)

target_link_libraries(finding_generators
        PUBLIC
        easy3d::core
        easy3d::renderer
        easy3d::viewer
        easy3d::gui
)

set(SPDLOG_DIR "${CMAKE_SOURCE_DIR}/spdlog")
message(STATUS "Using spdlog found at: ${SPDLOG_DIR}")
add_subdirectory(${SPDLOG_DIR})

target_link_libraries(finding_generators
        PUBLIC
        spdlog
)

set(JSON_DIR "${CMAKE_SOURCE_DIR}/json")
message(STATUS "Using json found at: ${JSON_DIR}")
add_subdirectory(${JSON_DIR})

target_link_libraries(finding_generators
        PUBLIC
        nlohmann_json::nlohmann_json
)

set(INIT_PY "${CMAKE_SOURCE_DIR}/package/__init__.py")
if(NOT EXISTS ${INIT_PY})
    file(WRITE ${INIT_PY} "")
endif()

add_custom_command(
        TARGET finding_generators
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_SOURCE_DIR}/package
        COMMAND ${CMAKE_COMMAND} -E copy
        $<TARGET_FILE:finding_generators>
        ${CMAKE_SOURCE_DIR}/package/
)