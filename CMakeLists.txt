cmake_minimum_required(VERSION 3.14)

message("${INCLUDE_DIRECTORIES}")

set(CMAKE_C_FLAGS "-std=c11 -fPIC ${CMAKE_C_FLAGS}")

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")

include_directories("/home/mikhail/Projects/density_fitting/libmol2_local_build/include")
link_directories(/home/mikhail/Projects/density_fitting/libmol2_local_build/lib)

add_executable(minimize minimize.c)
target_link_libraries(minimize mol2 jansson m)