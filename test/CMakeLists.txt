link_directories(${CHECK_LIBRARY_DIRS})

set(tests_list
        test_all)

# Configure data files
file(GLOB test_files "${CMAKE_CURRENT_SOURCE_DIR}/data/*")
foreach(filepath ${test_files})
    get_filename_component(filename ${filepath} NAME)
    configure_file(${filepath} ${CMAKE_CURRENT_BINARY_DIR}/${filename} COPYONLY)
endforeach()

set(CMAKE_EXE_LINKER_FLAGS "${CHECK_LDFLAGS_OTHER} ${CMAKE_EXE_LINKER_FLAGS}")

# Add and register tests
foreach(test ${tests_list})
    add_executable(${test}
            ${test}.c
            ${SRC_LIST})

    target_link_libraries(${test} PRIVATE
            ${CHECK_LIBRARIES}
            ${LIBRARY_LIST})

    target_include_directories(${test} PRIVATE
            ${PROJECT_SOURCE_DIR}/src
            ${CHECK_INCLUDE_DIRS}
            ${CHECK_INCLUDE_DIR}
            ${INCLUDE_LIST})

    add_test(${test} ${CMAKE_CURRENT_BINARY_DIR}/${test})
endforeach()
