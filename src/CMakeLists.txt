project(DIAAlignment)
include(CTest)
cmake_minimum_required(VERSION 2.8)
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

add_library(DIAAlignment SHARED ${SOURCE_FILES})

# Test executables
add_executable(FullTest test_alignment.cpp)
add_executable(runTest1 test_simpleFcn.cpp)
target_link_libraries(FullTest DIAAlignment)
target_link_libraries(runTest1 DIAAlignment)
add_test(RunTest FullTest)
add_test(RunTest1 runTest1)
