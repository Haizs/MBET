cmake_minimum_required(VERSION 3.16)
project(biclique)

set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -g -O3")

IF (CMAKE_BUILD_TYPE MATCHES Debug)
    add_definitions(-DDEBUG)
ENDIF ()

include_directories(./)

add_executable(biclique main.cpp)

add_executable(MBET main.cpp)
target_compile_definitions(MBET PRIVATE _MBET)

add_executable(MBETC main.cpp)
target_compile_definitions(MBETC PRIVATE _MBET ONLY_COUNTING)

add_executable(MBETR main.cpp)
target_compile_definitions(MBETR PRIVATE _MBET OUTPUT_RESULT)

add_executable(MBETM main.cpp)
target_compile_definitions(MBETM PRIVATE _MBET MBETM)

add_executable(MBETMC main.cpp)
target_compile_definitions(MBETMC PRIVATE _MBET MBETM ONLY_COUNTING)

add_executable(MBETMR main.cpp)
target_compile_definitions(MBETMR PRIVATE _MBET MBETM OUTPUT_RESULT)
