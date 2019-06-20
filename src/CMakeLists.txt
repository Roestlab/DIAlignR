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
CInterface.cpp
utils.cpp
)

add_library(DIAlign SHARED ${SOURCE_FILES})

# Test executables
add_executable(FullTest 
  test_alignment.cpp
  ${SOURCE_FILES})
add_test(RunTest FullTest) 

