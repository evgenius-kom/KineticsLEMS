cmake_minimum_required(VERSION 3.4)
set(CMAKE_CXX_STANDARD 17)

project(KineticsLEMS)

add_compile_options(-lzip -Wextra -Wall -Werror -pedantic-errors)

add_executable(${PROJECT_NAME} main.cpp)