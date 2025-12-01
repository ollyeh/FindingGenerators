cmake_minimum_required(VERSION 3.31)
project(FindingGenerators)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_OSX_DEPLOYMENT_TARGET 14.0)
include(FetchContent)

set(PYBIND11_VERSION 3.0.1)
FetchContent_Declare(
        pybind11
        GIT_REPOSITORY https://github.com/pybind/pybind11.git
        GIT_TAG master
        FIND_PACKAGE_ARGS ${PYBIND11_VERSION}
)
FetchContent_MakeAvailable(pybind11)

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

pybind11_add_module(finding_generators SHARED cpp_src/bindings.cpp)

execute_process(
        COMMAND ${Python3_EXECUTABLE} -m pybind11 --cmakedir
        OUTPUT_VARIABLE PYBIND11_CMAKE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
)

if (GRAPHICS)
    message(STATUS "Graphics enabled — building with Easy3D")
	target_compile_definitions(finding_generators PUBLIC GRAPHICS=1)
    set(EASY3D_VERSION 2.6.1)
    FetchContent_Declare(
            easy3d
            GIT_REPOSITORY https://github.com/LiangliangNan/Easy3D.git
            GIT_TAG main
            FIND_PACKAGE_ARGS ${EASY3D_VERSION}
    )
    set(EASY3D_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
    set(EASY3D_BUILD_TESTS    OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(easy3d)
    target_link_libraries(finding_generators
	    PRIVATE
            easy3d::core
            easy3d::renderer
            easy3d::viewer
            easy3d::gui
    )
    add_custom_command(
            TARGET finding_generators
            POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:easy3d::core>      ${CMAKE_SOURCE_DIR}
            COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:easy3d::renderer>  ${CMAKE_SOURCE_DIR}
            COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:easy3d::viewer>    ${CMAKE_SOURCE_DIR}
            COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:easy3d::gui>       ${CMAKE_SOURCE_DIR}
    )
else()
    message(STATUS "Graphics disabled — building WITHOUT Easy3D")
    target_compile_definitions(finding_generators PRIVATE GRAPHICS=0)
endif()

set(SPDLOG_VERSION 1.16.0)
FetchContent_Declare(
        spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG v1.x
        FIND_PACKAGE_ARGS ${SPDLOG_VERSION}
)
FetchContent_MakeAvailable(spdlog)

target_link_libraries(finding_generators
        PUBLIC
        spdlog
)

set(JSON_VERSION 3.12.0)
FetchContent_Declare(
        json
        GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_TAG develop
        FIND_PACKAGE_ARGS ${JSON_VERSION}
)
FetchContent_MakeAvailable(json)

target_link_libraries(finding_generators
        PUBLIC
        nlohmann_json::nlohmann_json
)

add_custom_command(
        TARGET finding_generators
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_SOURCE_DIR}
        COMMAND ${CMAKE_COMMAND} -E copy
        $<TARGET_FILE:finding_generators>
        ${CMAKE_SOURCE_DIR}
)
