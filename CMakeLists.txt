cmake_minimum_required(VERSION 3.9)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

cmake_policy(SET CMP0069 NEW)

project(NRGmin VERSION 1.0
        DESCRIPTION "Energy minimization utility"
        LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

option(BUILD_TESTS    "Build tests"                        ON)
option(USE_LTO        "Link Time Optimization"             ON)
option(USE_SANITIZER  "Build with address sanitizer"      OFF)
option(OPENMP         "Add OPENMP executable"             OFF)
option(NOE            "Build with NOE feature"             ON)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

if (UNIX AND CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-pg -Wall -Wshadow -Wpointer-arith -Wcast-qual -Winline -Werror -Wextra -Wfatal-errors -Wstrict-prototypes)
endif()

if (USE_SANITIZER)
    set(SANITIZER_FLAGS "-fno-omit-frame-pointer -fsanitize=address,undefined")
    add_compile_options(${SANITIZER_FLAGS})
    add_link_options(${SANITIZER_FLAGS})
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

# Read help text from the file and plug it in parse_options.c.in
file(STRINGS ${PROJECT_SOURCE_DIR}/src/usage.txt USAGE_LINES)
# Join lines using "\" character. Can't use string(JOIN ..) here,
# because it exists only in the newer versions of cmake
set(USAGE_TEXT "")
list(LENGTH USAGE_LINES USAGE_N_LINES)
math(EXPR USAGE_LAST_LINE "${USAGE_N_LINES} - 1")
foreach(LINE_ID RANGE ${USAGE_LAST_LINE})
    list(GET USAGE_LINES ${LINE_ID} LINE)
    set(USAGE_TEXT "${USAGE_TEXT}\\n${LINE}")
endforeach()
configure_file (
        "${PROJECT_SOURCE_DIR}/src/parse_options.c.in"
        "${PROJECT_BINARY_DIR}/parse_options.c"
)

set(SRC_LIST
        ${PROJECT_SOURCE_DIR}/src/setup.c
        ${PROJECT_SOURCE_DIR}/src/utils.c
        ${PROJECT_SOURCE_DIR}/src/energy.c
        ${PROJECT_BINARY_DIR}/parse_options.c)

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
