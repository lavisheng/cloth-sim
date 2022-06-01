cmake_minimum_required(VERSION 3.22)

project(cloth-sim)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#libigl setup

option(LIBIGL_WITH_OPENGL               "Use OpenGL"            ON)
option(LIBIGL_WITH_OPENGL_GLFW          "Use GLFW"              ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI    "Use ImGui"             ON)

include (libigl)

include_directories(${PROJECT_SOURCE_DIR}/include)

#include code
file(GLOB H_FILES ${PROJECT_SOURCE_DIR}/include/*.h)
file(GLOB CPP_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)

#exec building
add_executable(${PROJECT_NAME} main.cpp ${H_FILES} ${CPP_FILES})
target_link_libraries(${PROJECT_NAME} igl::core igl::opengl igl::opengl_glfw igl::opengl_glfw_imgui)
