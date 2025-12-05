cmake_minimum_required(VERSION 3.28)
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

pybind11_add_module(finding_generators_binding SHARED cpp_src/bindings.cpp)

if(APPLE)
    set_target_properties(finding_generators_binding PROPERTIES
	INSTALL_RPATH "@loader_path/easy3d"
        BUILD_RPATH "@loader_path/easy3d")
elseif(UNIX)
    set_target_properties(finding_generators_binding PROPERTIES
	INSTALL_RPATH "$ORIGIN/easy3d"
	BUILD_RPATH "$ORIGIN/easy3d")
endif()

execute_process(
        COMMAND ${Python3_EXECUTABLE} -m pybind11 --cmakedir
        OUTPUT_VARIABLE PYBIND11_CMAKE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
)

add_custom_command(
        TARGET finding_generators_binding
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_SOURCE_DIR}
        COMMAND ${CMAKE_COMMAND} -E copy
        $<TARGET_FILE:finding_generators_binding>
        ${CMAKE_SOURCE_DIR}
)

if (GRAPHICS)
	message(STATUS "Graphics enabled — building with Easy3D")
	target_compile_definitions(finding_generators_binding PUBLIC GRAPHICS=1)
	set(Easy3D_BUILD_TUTORIALS OFF CACHE BOOL "" FORCE)
	set(Easy3D_BUILD_TESTS OFF CACHE BOOL "" FORCE)
	set(Easy3D_BUILD_DOCUMENTATION OFF CACHE BOOL "" FORCE)

	set(EASY3D_VERSION 2.6.1)
	FetchContent_Declare(
		easy3d
        	GIT_REPOSITORY https://github.com/LiangliangNan/Easy3D.git
        	GIT_TAG main
        	FIND_PACKAGE_ARGS ${EASY3D_VERSION}
    		)
    	FetchContent_MakeAvailable(easy3d)
    	set(Easy3D_BUILD_TUTORIALS OFF CACHE BOOL "Disable tutorials" FORCE)
	set(Easy3D_BUILD_TESTS OFF CACHE BOOL "Disable tests" FORCE)
	add_custom_command(
		TARGET finding_generators_binding
		PRE_BUILD
		COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_SOURCE_DIR}/easy3d
		COMMAND ${CMAKE_COMMAND} -E copy_directory ${Easy3D_SOURCE_DIR}/resources ${CMAKE_SOURCE_DIR}/easy3d/resources
		)
	target_link_libraries(finding_generators_binding
		PRIVATE
        	easy3d::core
            	easy3d::renderer
            	easy3d::viewer
            	easy3d::gui
		easy3d::util
    		)
	set(NAMELIST "core;renderer;viewer;gui;fileio;algo;kdtree;util")
	if (APPLE)
		message(STATUS " - Symlinks for APPLE ")
	foreach(name IN LISTS NAMELIST)
		add_custom_command(
			TARGET finding_generators_binding
	    		POST_BUILD
	    		COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:easy3d::${name}> ${CMAKE_SOURCE_DIR}/easy3d
			COMMAND ${CMAKE_COMMAND} -E create_symlink $<TARGET_FILE_NAME:easy3d::${name}>    ${CMAKE_SOURCE_DIR}/easy3d/libeasy3d_${name}.2.dylib
			)
	endforeach()
	elseif(UNIX)
		message(STATUS " - Symlinks for UNIX ")
	else()
		message(STATUS " - Symlinks for else ")
	endif()
else()
    message(STATUS "Graphics disabled — building WITHOUT Easy3D")
    target_compile_definitions(finding_generators_binding PRIVATE GRAPHICS=0)
endif()

set(SPDLOG_VERSION 1.16.0)
FetchContent_Declare(
        spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG v1.x
        FIND_PACKAGE_ARGS ${SPDLOG_VERSION}
)
FetchContent_MakeAvailable(spdlog)
set_property(TARGET spdlog PROPERTY POSITION_INDEPENDENT_CODE ON)
target_link_libraries(finding_generators_binding
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

target_link_libraries(finding_generators_binding
        PUBLIC
        nlohmann_json::nlohmann_json
)

