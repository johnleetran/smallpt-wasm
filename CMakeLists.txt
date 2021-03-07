cmake_minimum_required(VERSION 3.10)

#set project name
project(smallptWASM VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)


include_directories("./include")
add_executable(smallptWASM smallpt-dynamic.cpp)
