project(DIAlign)
include(CTest)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)

aux_source_directory(. SRC_LIST)
set (CMAKE_CXX_STANDARD 11)
set(SOURCE_FILES 
affinealignment.cpp
affinealignobj.cpp
alignment.cpp
chromSimMatrix.cpp
constrainMat.cpp
gapPenalty.cpp
utils.cpp
)

add_library(DIAlign SHARED ${SOURCE_FILES})

# Test executables
add_executable(FullTest test_alignment.cpp)
add_executable(runAlignment run_alignment.cpp)
target_link_libraries(FullTest DIAlign)
target_link_libraries(runAlignment DIAlign)
add_test(RunTest FullTest) 

