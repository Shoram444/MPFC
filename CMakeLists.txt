cmake_minimum_required(VERSION 2.8.9)
project(MPFC)

find_package(ROOT)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(include ${ROOT_INCLUDE_DIRS})
file(GLOB SOURCES "src/*.cpp" "src/dicts/*.cpp")


add_library(MPFC SHARED ${SOURCES} )
target_link_libraries(MPFC PUBLIC ${ROOT_LIBRARIES})






 
