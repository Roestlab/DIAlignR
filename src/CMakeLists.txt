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
simpleFcn.cpp
)

add_library(DIAAlignment SHARED ${SOURCE_FILES})
# SHARED libraries are linked dynamically and loaded at runtime. Other options are
# STATIC or MODULE

# Test executables
add_executable(FullTest test_alignment.cpp)
add_executable(runTest1 test_simpleFcn.cpp)
add_executable(runTest2 test_utils.cpp)
target_link_libraries(FullTest DIAAlignment)
target_link_libraries(runTest1 DIAAlignment)
target_link_libraries(runTest2 DIAAlignment)
add_test(RunTest FullTest)
add_test(RunTest1 runTest1)
add_test(RunTest2 runTest2)

# https://stackoverflow.com/questions/39598323/how-to-properly-link-libraries-with-cmake
