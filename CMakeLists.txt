cmake_minimum_required(VERSION 3.10)

# Set the project name
project(CMS-BWT)

add_subdirectory(libsais)

# Add the executable
add_executable(cmsbwt CMS-BWT-functions.cpp CMS-BWT.h match.h main.cpp)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native")

target_include_directories(cmsbwt PUBLIC libsais/include)
target_link_libraries(cmsbwt libsais)