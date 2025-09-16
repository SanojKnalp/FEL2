#Cmake file to setup the project

#Cmake minimum required version
cmake_minimum_required(VERSION 4.0.2)

#Project name and languages we are using. Currently C++ and C. We need C for PaStiX later.
project(FEL2 LANGUAGES CXX C)

#Print CMake version to screen
message(Status "Using CMake ${CMAKE_VERSION}$")

#Set C++ Standards
set(CMAKE_CXX_STANDARD 20)
#Enforce the standard
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

#We will use a debuggable version cause we are working on the project
set(DEFAULT_BUILD_TYPE "Debug")

#search for all cc files in 
file(GLOB_RECURSE SOURCES "src/*.cc")

#generate executable
add_executable(FEL2 ${SOURCES})

target_include_directories(FEL2 PRIVATE src)
