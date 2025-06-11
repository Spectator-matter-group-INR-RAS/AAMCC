cmake_minimum_required(VERSION 3.16)
project(FermiBreakUp)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wpedantic -Werror")

find_package(CLHEP REQUIRED)

set(Utilities
        Utilities/DataTypes.cpp
        Utilities/IntegerPartition.cpp
        Utilities/NucleiProperties/VNucleiProperties.cpp
        Utilities/NucleiProperties/FermiNucleiProperties.cpp
        Utilities/NucleiProperties/Builder/CSVBuilder.cpp
        Utilities/NucleiProperties/Builder/DefaultBuilder.cpp
        Utilities/ConfigurationProperties.cpp
        Utilities/Randomizer.cpp)

set(PhaseSpace
        PhaseDecay/VDecay.cpp
        PhaseDecay/Decay.cpp
        PhaseDecay/KopylovDecay.cpp)

set(FragmentPool
        FragmentPool/FermiFragmentPool.cpp
        FragmentPool/Builder/DefaultBuilder.cpp
        FragmentPool/Fragments/FermiFragment.cpp
        FragmentPool/Fragments/StableFermiFragment.cpp
        FragmentPool/Fragments/UnstableFermiFragment.cpp
        FragmentPool/Fragments/B9FermiFragment.cpp
        FragmentPool/Fragments/Be8FermiFragment.cpp
        FragmentPool/Fragments/Li5FermiFragment.cpp
        FragmentPool/Fragments/He5FermiFragment.cpp)

set(FermiBreakUpModel
        VFermiBreakUp.cpp
        FermiParticle.cpp
        FermiSplit.cpp
        FermiConfigurations.cpp
        CachedFermiConfigurations.cpp
        FermiBreakUp.cpp)

set(Sources ${FermiBreakUpModel} ${FragmentPool} ${Utilities} ${PhaseSpace})

add_library(FermiBreakUp STATIC ${Sources})

target_link_libraries(FermiBreakUp CLHEP::CLHEP)

target_include_directories(FermiBreakUp PUBLIC ${PROJECT_SOURCE_DIR})
