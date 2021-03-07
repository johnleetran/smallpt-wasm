cmake_minimum_required(VERSION 3.10)

#set project name
project(smallptwasm VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)


include_directories("./include")
# add_executable(smallptwasm smallpt-dynamic.cpp)
add_library( smallptwasm STATIC smallpt-dynamic.cpp )