cmake_minimum_required(VERSION 3.16)
project(NuFitSampler CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(ROOT REQUIRED COMPONENTS Core RIO Tree)
find_package(Threads REQUIRED)

add_library(nufitcore
  src/GridTable.cc
  src/NuFitModel.cc
  src/DirectSampler.cc
  src/PTSampler.cc
  src/NuFitTarget.cc
  src/Output.cc)
target_include_directories(nufitcore PUBLIC include)
target_link_libraries(nufitcore PUBLIC ROOT::Core ROOT::RIO ROOT::Tree Threads::Threads)
target_compile_options(nufitcore PRIVATE -Wall -Wextra)

foreach(app sample_direct sample_pt)
  add_executable(${app} apps/${app}.cc)
  target_link_libraries(${app} PRIVATE nufitcore)
endforeach()

enable_testing()
add_executable(test_core tests/test_core.cc)
target_link_libraries(test_core PRIVATE nufitcore)
add_test(NAME core COMMAND test_core)
add_executable(test_nufit tests/test_nufit.cc)
target_link_libraries(test_nufit PRIVATE nufitcore)
add_test(NAME nufit COMMAND test_nufit ${CMAKE_SOURCE_DIR}/data/nufit)
set_tests_properties(nufit PROPERTIES SKIP_RETURN_CODE 77)
