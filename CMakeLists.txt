cmake_minimum_required(VERSION 3.9)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

if(CMAKE_VERSION VERSION_GREATER 3.9)
    cmake_policy(SET CMP0069 NEW)
endif()

project(NRGmin VERSION 1.0
        DESCRIPTION "Energy minimization utility"
        LANGUAGES C)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(BUILD_TESTS    "Build tests"                        ON)
option(USE_LTO        "Link Time Optimization"             ON)
option(USE_SANITIZER  "Build with address sanitizer"      OFF)
option(OPENMP         "Add OPENMP executable"             OFF)
option(NOE            "Build with NOE feature"             ON)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

if (UNIX)
    add_compile_options(-std=c11)
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        add_compile_options(-pg -Wall -Wshadow -Wpointer-arith -Wcast-qual -Winline -Werror -Wextra -Wfatal-errors -Wstrict-prototypes)
    endif()
endif()

if (USE_SANITIZER)
    set(SANITIZER_FLAGS "-fno-omit-frame-pointer -fsanitize=address,undefined")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${SANITIZER_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}  ${SANITIZER_FLAGS}")
endif()

if (NOE)
    add_definitions(-D NOE)
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

if (USE_LTO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT result)
    if(NOT result)
        set(USE_LTO OFF)
    endif()
endif()

if (USE_LTO)
    set_target_properties(nrgmin PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

target_include_directories(nrgmin
        PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>"
        "$<INSTALL_INTERFACE:${INCLUDE_INSTALL_DIR}>"
        PRIVATE ${INCLUDE_LIST})

target_link_libraries(nrgmin ${LIBRARY_LIST})

if (OPENMP)
    find_package(OpenMP REQUIRED)

    add_executable(nrgmin.omp ${PROJECT_SOURCE_DIR}/src/minimize.c ${SRC_LIST})

    if (USE_LTO)
        set_target_properties(nrgmin.omp PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
    endif()

    target_include_directories(nrgmin.omp
        PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>"
        "$<INSTALL_INTERFACE:${INCLUDE_INSTALL_DIR}>"
        PRIVATE ${INCLUDE_LIST})

    target_link_libraries(nrgmin.omp ${LIBRARY_LIST} OpenMP::OpenMP_C)

    target_compile_definitions(nrgmin.omp PUBLIC OPENMP)
endif()

if (BUILD_TESTS)
    find_package(Check REQUIRED)
    enable_testing()
    add_subdirectory("test")
endif()