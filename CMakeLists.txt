#cmake -DCMAKE_CXX_FLAGS=-pg -DCMAKE_EXE_LINKER_FLAGS=-pg -DCMAKE_SHARED_LINKER_FLAGS=-pg -DCMAKE_BUILD_TYPE=Debug ../..
#cmake -DCMAKE_BUILD_TYPE=Release ../..
cmake_minimum_required(VERSION 3.10.0)
project(pSims-SphericalWave VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(exec "${PROJECT_SOURCE_DIR}/src/main.cpp"
                    "${PROJECT_SOURCE_DIR}/src/particle.cpp"
                    "${PROJECT_SOURCE_DIR}/src/wave.cpp"
                    "${PROJECT_SOURCE_DIR}/src/simulation.cpp")

set_target_properties(exec PROPERTIES DEBUG_POSTFIX "${CMAKE_DEBUG_POSTFIX}")

target_include_directories(exec PUBLIC 
            "${PROJECT_SOURCE_DIR}/headers/"
)

find_package(OpenMP COMPONENTS CXX REQUIRED)
find_package(HDF5 COMPONENTS CXX REQUIRED)
target_link_libraries(exec 
            OpenMP::OpenMP_CXX 
            HDF5::HDF5
)