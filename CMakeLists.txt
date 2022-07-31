cmake_minimum_required(VERSION 3.0)
project(PIPG VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CXX_FLAGS "-O3")
SET(CMAKE_CXX_FLAGS "-Wall")

#add_library(PIPG SHARED)
add_executable(PIPG src/test.cpp src/pipg.cpp)
target_include_directories(PIPG PUBLIC include)
#target_sources(PIPG PRIVATE src/pipg.cpp)

target_include_directories(PIPG PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(PIPG PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/lib/eigen)