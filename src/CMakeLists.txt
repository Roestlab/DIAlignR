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
add_executable(runTest3 test_gapPenalty.cpp)
add_executable(runTest4 test_constrainMat.cpp)
add_executable(runTest5 test_chromSimMatrix.cpp)

target_link_libraries(FullTest DIAAlignment)
target_link_libraries(runTest1 DIAAlignment)
target_link_libraries(runTest2 DIAAlignment)
target_link_libraries(runTest3 DIAAlignment)
target_link_libraries(runTest4 DIAAlignment)
target_link_libraries(runTest5 DIAAlignment)

add_test(RunTest FullTest)
add_test(RunTest1 runTest1)
add_test(RunTest2 runTest2)
add_test(RunTest3 runTest3)
add_test(RunTest4 runTest4)
add_test(RunTest5 runTest5)

# https://stackoverflow.com/questions/39598323/how-to-properly-link-libraries-with-cmake
