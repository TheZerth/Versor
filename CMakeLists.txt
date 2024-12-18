cmake_minimum_required(VERSION 3.29.6)
project(Versor)

set(CMAKE_CXX_STANDARD 23)

add_executable(Versor main.cpp
        libs/Versor.cpp
        libs/Versor.h)
