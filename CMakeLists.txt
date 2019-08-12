cmake_minimum_required(VERSION 3.5)
project(energymin LANGUAGES C)
set(energymin_version_major 0)
set(energymin_version_minor 0)
set(energymin_version "${energymin_version_major}.${energymin_version_minor}")
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

set(CMAKE_C_FLAGS "-std=c11 -fPIC ${CMAKE_C_FLAGS}")
set(CMAKE_C_FLAGS_DEBUG "-Wall -pg ${CMAKE_C_FLAGS_DEBUG}")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")

include_directories($ENV{HOME}/include)
link_directories($ENV{HOME}/lib)

find_package(GSL REQUIRED)

add_executable(minimize minimize.c)
target_link_libraries(minimize mol2 nmrgrad ${GSL_LIBRARIES} jansson m)

find_package(Check)
if(CHECK_FOUND)
    enable_testing()
    add_subdirectory(tests)
else()
    message("Could not find Check for testing")
endif()