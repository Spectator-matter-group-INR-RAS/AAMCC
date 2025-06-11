cmake_minimum_required(VERSION 3.16)
project(AAMCCFermiBreakUp)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wpedantic -Werror")

find_package(Geant4 REQUIRED)

include(${Geant4_USE_FILE})

add_library(AAMCCFermiBreakUp STATIC AAMCCFermiBreakUp.cc)

add_subdirectory(lib)
target_link_libraries(AAMCCFermiBreakUp FermiBreakUp ${Geant4_LIBRARIES})

