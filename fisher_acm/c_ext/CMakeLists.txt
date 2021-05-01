cmake_minimum_required(VERSION 3.19)
project(fexact C)

set(CMAKE_C_STANDARD 11)

add_executable(fexact main.c fexact.c helpers.c helpers.h fexact.h)