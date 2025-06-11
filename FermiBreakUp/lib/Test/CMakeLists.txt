cmake_minimum_required(VERSION 3.16)
project(FermiBreakUp)

set(CMAKE_CXX_STANDARD 20)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wpedantic -Werror")

set(LibPath ${CMAKE_CURRENT_SOURCE_DIR}/..)

find_package(CLHEP REQUIRED)

set(Tables
        ${LibPath}/TableValues/NucleiPropertiesTable.cpp
        ${LibPath}/TableValues/NucleiPropertiesTableAME12.cpp
)

set(Utilities
        ${LibPath}/Utilities/DataTypes.cpp
        ${LibPath}/Utilities/IntegerPartition.cpp
        ${LibPath}/Utilities/NucleiProperties/VNucleiProperties.cpp
        ${LibPath}/Utilities/NucleiProperties/Impl/FermiNucleiProperties.cpp
        ${LibPath}/Utilities/NucleiProperties/Impl/FastNucleiProperties.cpp
        ${LibPath}/Utilities/NucleiProperties/Impl/NativeNucleiProperties.cpp
        ${LibPath}/Utilities/NucleiProperties/DataStorage/CSVNuclearMass.cpp
        ${LibPath}/Utilities/NucleiProperties/DataStorage/DefaultNuclearMass.cpp
        ${LibPath}/Utilities/ConfigurationProperties.cpp
        ${LibPath}/Utilities/Randomizer.cpp
)

set(PhaseSpace
        ${LibPath}/PhaseDecay/VDecay.cpp
        ${LibPath}/PhaseDecay/Decay.cpp
        ${LibPath}/PhaseDecay/KopylovDecay.cpp
)

set(FragmentPool
        ${LibPath}/FragmentPool/FragmentsStorage.cpp
        ${LibPath}/FragmentPool/DataSource/DefaultPoolSource.cpp
        ${LibPath}/FragmentPool/Fragments/FermiFragment.cpp
        ${LibPath}/FragmentPool/Fragments/StableFermiFragment.cpp
        ${LibPath}/FragmentPool/Fragments/UnstableFermiFragment.cpp
        ${LibPath}/FragmentPool/Fragments/B9FermiFragment.cpp
        ${LibPath}/FragmentPool/Fragments/Be8FermiFragment.cpp
        ${LibPath}/FragmentPool/Fragments/Li5FermiFragment.cpp
        ${LibPath}/FragmentPool/Fragments/He5FermiFragment.cpp
)

set(Configurations
        ${LibPath}/Configurations/FermiConfigurations.cpp
        ${LibPath}/Configurations/CachedFermiConfigurations.cpp
        ${LibPath}/Configurations/FastFermiConfigurations.cpp
)

set(FermiBreakUpModel
        ${LibPath}/VFermiBreakUp.cpp
        ${LibPath}/FermiParticle.cpp
        ${LibPath}/FermiSplit.cpp
        ${LibPath}/FermiBreakUp.cpp
        ${LibPath}/VFermiConfigurations.cpp
)

set(Sources ${FermiBreakUpModel} ${FragmentPool} ${Utilities} ${PhaseSpace} ${Tables} ${Configurations})

enable_testing()
find_package(GTest REQUIRED)
include(GoogleTest)

set(Tests
        PartitionTest.cpp
        TableTests.cpp
        BreakUpTest.cpp
)

add_executable(FermiTest ${Tests} ${Sources})

target_link_libraries(FermiTest CLHEP::CLHEP)
target_link_libraries(FermiTest GTest::GTest GTest::Main)

target_include_directories(FermiTest PUBLIC ${LibPath})

gtest_discover_tests(FermiTest)
