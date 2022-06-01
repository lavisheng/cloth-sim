cmake_minimum_required(VERSION 3.22)

project(cloth-sim)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#libigl setup

option(LIBIGL_OPENGL               "Use OpenGL"            ON)
option(LIBIGL_GLFW                 "Use glfw"              ON)
find_package(LIBIGL)
#include (libigl)

include_directories(${PROJECT_SOURCE_DIR}/include)

#include code
file(GLOB H_FILES ${PROJECT_SOURCE_DIR}/include/*.h)
file(GLOB CPP_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)

#exec building
add_executable(${PROJECT_NAME} main.cpp ${H_FILES} ${CPP_FILES})
target_link_libraries(${PROJECT_NAME} igl::core igl::glfw igl::opengl)
