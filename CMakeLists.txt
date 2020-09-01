cmake_minimum_required(VERSION 2.8.9)
project(MPFC)


include_directories(include /home/shoram/Work/ROOT/install/include)
file(GLOB SOURCES "src/*.cpp")

add_library(MPFC SHARED ${SOURCES} ${ROOT_Physics_LIBRARY})







 
