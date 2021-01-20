project(DIAlign)
include(CTest)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
if (MSVC)
  # MSVC specific flags
else()
  SET(CMAKE_CXX_FLAGS  "-march=native")  # about a 10% improvement
endif()

aux_source_directory(. SRC_LIST)
set (CMAKE_CXX_STANDARD 11)
set(SOURCE_FILES
src/affinealignment.cpp
src/affinealignobj.cpp
src/alignment.cpp
src/chromSimMatrix.cpp
src/constrainMat.cpp
src/gapPenalty.cpp
src/utils.cpp
src/simpleFcn.cpp
src/integrateArea.cpp
src/PeakIntegrator.cpp
src/MSChromatogram.cpp
src/ChromatogramPeak.cpp
src/DPosition.cpp
src/spline.cpp
src/miscell.cpp
)

add_library(DIAAlignment ${SOURCE_FILES})
target_compile_definitions(DIAAlignment PRIVATE -DDIALIGN_PURE_CPP=On)
# SHARED libraries are linked dynamically and loaded at runtime. Other options are
# STATIC or MODULE

# Actual executables
add_executable(runAlignment src/run_alignment.cpp)
target_link_libraries(runAlignment DIAAlignment)
target_compile_definitions(runAlignment PRIVATE -DDIALIGN_PURE_CPP=On)

# Test executables
add_executable(runTest1 src/test/test_simpleFcn.cpp)
add_executable(runTest2 src/test/test_utils.cpp)
add_executable(runTest3 src/test/test_gapPenalty.cpp)
add_executable(runTest4 src/test/test_constrainMat.cpp)
add_executable(runTest5 src/test/test_chromSimMatrix.cpp)
add_executable(runTest6 src/test/test_alignment.cpp)
add_executable(runTest7 src/test/test_affinealignobj.cpp)
add_executable(runTest8 src/test/test_affinealignment.cpp)
add_executable(runTest9 src/test/test_integrateArea.cpp)
add_executable(runTest10 src/test/test_miscell.cpp)

set(LIST_TESTS
runTest1
runTest2
runTest3
runTest4
runTest5
runTest6
runTest7
runTest8
runTest9
runTest10
)

foreach(TEST ${LIST_TESTS})
  target_compile_definitions(${TEST} PRIVATE -DDIALIGN_PURE_CPP=On)
  target_link_libraries(${TEST} DIAAlignment)
  add_test(${TEST} ${TEST})
endforeach()

# https://stackoverflow.com/questions/39598323/how-to-properly-link-libraries-with-cmake
