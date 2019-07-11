cmake_minimum_required(VERSION 3.5)

message("${INCLUDE_DIRECTORIES}")

set(CMAKE_C_FLAGS "-std=c11 -fPIC ${CMAKE_C_FLAGS}")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")

include_directories(/home/mikhail/include)
link_directories(/home/mikhail/lib)

find_package(PkgConfig REQUIRED)
pkg_check_modules(GSL REQUIRED gsl)

add_executable(minimize minimize.c)
target_link_libraries(minimize mol2 nmrgrad ${GSL_LIBRARIES} jansson m)
