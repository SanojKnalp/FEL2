cmake_minimum_required(VERSION 3.31)
project(COME LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Armadillo (via vcpkg)
find_package(Armadillo CONFIG QUIET)

if(NOT Armadillo_FOUND)
    find_package(Armadillo REQUIRED)
endif()

# Core library
add_library(come_lib)

target_link_libraries(come_lib
    PUBLIC armadillo 
)


target_include_directories(come_lib
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/source
)

add_subdirectory(source)

# Main executable
add_executable(COME main.cpp)
target_link_libraries(COME PRIVATE come_lib)

# Tests
enable_testing()
add_subdirectory(UnitTests)
add_subdirectory(SystemTests)