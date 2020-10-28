cmake_minimum_required(VERSION 3.13)
project(Problem_komw)

set(CMAKE_CXX_STANDARD 17)

add_compile_options(-Wall -Wextra -Werror -Wpedantic -pedantic-errors)

include_directories(
        include
)

set(SOURCE_FILES
        src/tsp.cpp)


add_executable(Problem_komw main.cpp ${SOURCE_FILES})
