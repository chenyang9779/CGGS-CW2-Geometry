# Set the project name
cmake_minimum_required(VERSION 3.20)
project(Section2)

# Add the executable
add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME} common_includes polyscope)

