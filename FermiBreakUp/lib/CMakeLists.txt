cmake_minimum_required(VERSION 3.16)
project(FermiBreakUp)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wpedantic -Werror")

find_package(CLHEP REQUIRED)

set(Tables
        TableValues/NucleiPropertiesTable.cpp
        TableValues/NucleiPropertiesTableAME12.cpp
)

set(Utilities
        Utilities/DataTypes.cpp
        Utilities/IntegerPartition.cpp
        Utilities/NucleiProperties/VNucleiProperties.cpp
        Utilities/NucleiProperties/Impl/FermiNucleiProperties.cpp
        Utilities/NucleiProperties/Impl/FastNucleiProperties.cpp
        Utilities/NucleiProperties/Impl/NativeNucleiProperties.cpp
        Utilities/NucleiProperties/DataStorage/CSVNuclearMass.cpp
        Utilities/NucleiProperties/DataStorage/DefaultNuclearMass.cpp
        Utilities/ConfigurationProperties.cpp
        Utilities/Randomizer.cpp
)

set(PhaseSpace
        PhaseDecay/VDecay.cpp
        PhaseDecay/Decay.cpp
        PhaseDecay/KopylovDecay.cpp
)

set(FragmentPool
        FragmentPool/FragmentsStorage.cpp
        FragmentPool/DataSource/DefaultPoolSource.cpp
        FragmentPool/Fragments/FermiFragment.cpp
        FragmentPool/Fragments/StableFermiFragment.cpp
        FragmentPool/Fragments/UnstableFermiFragment.cpp
        FragmentPool/Fragments/B9FermiFragment.cpp
        FragmentPool/Fragments/Be8FermiFragment.cpp
        FragmentPool/Fragments/Li5FermiFragment.cpp
        FragmentPool/Fragments/He5FermiFragment.cpp
)

set(Configurations
        Configurations/FermiConfigurations.cpp
        Configurations/CachedFermiConfigurations.cpp
        Configurations/FastFermiConfigurations.cpp
)

set(FermiBreakUpModel
        VFermiBreakUp.cpp
        FermiParticle.cpp
        FermiSplit.cpp
        FermiBreakUp.cpp
        VFermiConfigurations.cpp
)

set(Sources ${FermiBreakUpModel} ${FragmentPool} ${Utilities} ${PhaseSpace} ${Tables} ${Configurations})

add_library(FermiBreakUp STATIC ${Sources})

target_link_libraries(FermiBreakUp CLHEP::CLHEP)

target_include_directories(FermiBreakUp PUBLIC ${PROJECT_SOURCE_DIR})
