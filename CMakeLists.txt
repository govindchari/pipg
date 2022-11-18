cmake_minimum_required(VERSION 3.0)
project(PIPG VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_FLAGS "-O3")
SET(CMAKE_CXX_FLAGS "-Wall")

add_library(PIPG SHARED)
target_include_directories(PIPG PUBLIC include)
target_sources(PIPG PRIVATE src/pipg.cpp)
target_include_directories(PIPG PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(PIPG PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/lib/eigen)

#Testing
if (ENABLE_TESTING)
    message("PIPG Testing Enabled")
    include_directories(tests)
    add_executable(run_tests tests/run_tests.cpp)
    target_include_directories(run_tests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_include_directories(run_tests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/lib/eigen)
    target_include_directories(run_tests PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/tests)

    if (ENABLE_COVERAGE)
        SET(CMAKE_CXX_FLAGS "--coverage -O0")
        target_link_libraries(run_tests PIPG --coverage)
    else()
        target_link_libraries(run_tests PIPG) 
    endif()
endif()
