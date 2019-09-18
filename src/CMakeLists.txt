project(DIAlign)
include(CTest)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
SET(CMAKE_CXX_FLAGS  "-march=native")  # about a 10% improvement

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
target_compile_definitions(DIAAlignment PRIVATE -DPURE_CPP=On)
# SHARED libraries are linked dynamically and loaded at runtime. Other options are
# STATIC or MODULE

# Actual executables
add_executable(runAlignment run_alignment.cpp)
target_link_libraries(runAlignment DIAAlignment)
target_compile_definitions(runAlignment PRIVATE -DPURE_CPP=On)

# Test executables
add_executable(runTest1 test_simpleFcn.cpp)
add_executable(runTest2 test_utils.cpp)
add_executable(runTest3 test_gapPenalty.cpp)
add_executable(runTest4 test_constrainMat.cpp)
add_executable(runTest5 test_chromSimMatrix.cpp)
add_executable(runTest6 test_alignment.cpp)
add_executable(runTest7 test_affinealignobj.cpp)
add_executable(runTest8 test_affinealignment.cpp)

set(LIST_TESTS
runTest1
runTest2
runTest3
runTest4
runTest5
runTest6
runTest7
runTest8
)

foreach(TEST ${LIST_TESTS})
  target_compile_definitions(${TEST} PRIVATE -DPURE_CPP=On)
  target_link_libraries(${TEST} DIAAlignment)
  add_test(${TEST} ${TEST})
endforeach()

# https://stackoverflow.com/questions/39598323/how-to-properly-link-libraries-with-cmake
