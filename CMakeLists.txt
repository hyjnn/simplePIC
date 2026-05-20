cmake_minimum_required(VERSION 3.28)
project(myPIC VERSION 0.1 DESCRIPTION "Practice PIC simulation code" LANGUAGES CXX)

add_library(MDVeclib INTERFACE)
target_include_directories(MDVeclib INTERFACE src/)
target_compile_features(MDVeclib INTERFACE cxx_std_23)

add_library(FDTDlib SHARED src/FDTD.cpp)
target_include_directories(FDTDlib PUBLIC src/)
target_compile_features(FDTDlib PUBLIC cxx_std_23)
target_link_libraries(FDTDlib PUBLIC MDVeclib)

add_executable(PIC src/main.cpp)
target_include_directories(PIC PRIVATE src/) # Redundant right now but kept here in case I move the library headers
target_compile_features(PIC PRIVATE cxx_std_23)
target_link_libraries(PIC PRIVATE FDTDlib)