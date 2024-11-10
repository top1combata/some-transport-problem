cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 20)

project(main)

add_subdirectory(submodules/quadpackcpp)
add_subdirectory(submodules/eigen)


add_executable(main
    src/main.cpp
    src/DualFunction.h
    src/ProblemFunctional.cpp
    src/ProblemFunctional.h
    src/Solver.cpp
    src/Solver.h
)

target_link_libraries(main quadpackcpp)
target_link_libraries(main eigen)
