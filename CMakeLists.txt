cmake_minimum_required(VERSION 3.13)
project(proj)

set(CMAKE_CXX_STANDARD 11)

add_executable(proj main.cpp Simulator.cpp Simulator.h dCmpUtil.h)

add_executable(test main2.cpp)