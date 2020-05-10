cmake_minimum_required(VERSION 3.1...3.15)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

project(EnergyMin VERSION 1.0
        DESCRIPTION "Energy minimization utility"
        LANGUAGES C)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(LTO "Link Time Optimization" ON)
option(BUILD_TESTS "Build tests" ON)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-Wall -Wshadow -Wpointer-arith -Wcast-qual -Winline -Werror -Wextra -Wfatal-errors -Wstrict-prototypes)
endif()

find_package(Jansson REQUIRED)

find_package(mol2 REQUIRED)

set(LIBRARY_LIST
        "${MOL2_LIBRARIES}"
        "${JANSSON_LIBRARIES}"
        m)

set(INCLUDE_LIST
        "${MOL2_INCLUDE_DIRS}"
        "${JANSSON_INCLUDE_DIRS}")

set(SRC_LIST
        ${PROJECT_SOURCE_DIR}/src/setup.c
        ${PROJECT_SOURCE_DIR}/src/parse_options.c
        ${PROJECT_SOURCE_DIR}/src/utils.c
        ${PROJECT_SOURCE_DIR}/src/energy.c)

add_executable(nrgmin ${PROJECT_SOURCE_DIR}/src/minimize.c ${SRC_LIST})

if (LTO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT result)
    if(result)
        set_target_properties(nrgmin PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif()
endif()

target_include_directories(nrgmin
        PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>"
        "$<INSTALL_INTERFACE:${INCLUDE_INSTALL_DIR}>"
        PRIVATE ${INCLUDE_LIST})

target_link_libraries(nrgmin ${LIBRARY_LIST})

if (BUILD_TESTS)
    find_package(Check REQUIRED)
    enable_testing()
    add_subdirectory("test")
endif()