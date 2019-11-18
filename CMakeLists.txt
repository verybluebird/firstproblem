cmake_minimum_required(VERSION 3.14)
project(RNA)

add_subdirectory(googletest)

include_directories(googletest/googletest/include)

include_directories(googletest/googlemock/include)

set(CMAKE_CXX_STANDARD 14)

add_executable(RNA main.cpp RNA.cpp RNA.h TEST.cpp)

target_link_libraries(RNA gtest gtest_main)

target_link_libraries(RNA gmock gmock_main)
