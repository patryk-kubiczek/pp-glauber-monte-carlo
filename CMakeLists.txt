cmake_minimum_required(VERSION 3.1)
project(PP_MC)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(LIBS -lm -lgsl -lgslcblas )

set(SOURCE_FILES main.cpp pp.cpp)
add_executable(PP_MC ${SOURCE_FILES})

target_link_libraries(PP_MC  ${LIBS})


