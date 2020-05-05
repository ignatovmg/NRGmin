cmake_minimum_required(VERSION 3.1...3.15)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

project(EnergyMin VERSION 0.0
        DESCRIPTION "Energy minimization tool"
        LANGUAGES C)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(LTO "Link Time Optimization" ON)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

#add_compile_options(-Wall -Wshadow -Wpointer-arith -Wcast-qual -Winline -Werror -Wextra -Wfatal-errors -Wstrict-prototypes)

include(CMakePrintHelpers)

cmake_print_variables(CMAKE_MODULE_PATH CMAKE_MODULE_PATH)

find_package(Jansson REQUIRED)

find_package(GSL REQUIRED)

find_package(mol2 REQUIRED)

add_executable(minimize src/minimize.c src/potentials.c src/potentials.h src/parse_options.c src/parse_options.h src/utils.h)

if (LTO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT result)
    if(result)
        set_target_properties(minimize PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif()
endif()

target_include_directories(minimize
        PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>"
        "$<INSTALL_INTERFACE:${INCLUDE_INSTALL_DIR}>"
        PRIVATE "${MOL2_INCLUDE_DIRS} ${JANSSON_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS} ${GSL_INCLUDE_DIR}")

target_link_libraries(minimize ${MOL2_LIBRARIES} ${GSL_LIBRARIES} ${JANSSON_LIBRARIES} m)
