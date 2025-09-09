cmake_minimum_required(VERSION 3.12.0)
project(singlePhaseLTORC)

# Set the installation prefix
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR} CACHE PATH "Installation directory" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

## Default to Release build if no build type is specified
## Check if CMAKE_BUILD_TYPE is set; if not, set it to Release
#if(NOT CMAKE_BUILD_TYPE)
#    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Build type (default: Debug)" FORCE)
#endif()
#message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")


# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# execute_process(COMMAND bash -c "source /opt/scientific/intel/bin/compilervars.sh intel64")

set(fmt_DIR /backup5/weiye/SLTORC/fmt/fmt-install/lib/cmake/fmt) # path/to/ .cmake file for FMT
#set(ENV{PKG_CONFIG_PATH} "/opt/scientific/cantera-2.4_gnu_blas/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}") # path/to/ .pc file for CANTERA

set(ENV{PKG_CONFIG_PATH} "/backup/weiye/cantera_2.4.0/installDir/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
# Find necessary packages
find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(fmt REQUIRED)

# Use pkg-config to find GSL and Cantera
find_package(PkgConfig REQUIRED)
pkg_check_modules(GSL REQUIRED IMPORTED_TARGET gsl)
pkg_check_modules(CANTERA REQUIRED IMPORTED_TARGET cantera)


# Define path to header files and libraries
#set(IDA_DIR /home/weiye/src/sundials-5.7.0/installDir/lib)
set(IDA_DIR /home/weiye/src/sundials-3.1.1/installDir/lib)
set(COOLPROP_DIR /backup/weiye/CoolProp/install_root/shared_library/Linux)
set(FMT_DIR /backup/weiye/fmt/lib)

set(IDA_INCLUDES ${IDA_DIR}/../include)
set(COOLPROP_INCLUDES ${COOLPROP_DIR}/../)
set(FMT_INCLUDES ${FMT_DIR}/../include)

#find_library(SUNDIALS_NVECOPENMP_LIB NAMES sundials_nvecopenmp PATHS ${IDA_DIR} NO_DEFAULT_PATH)
#find_library(SUNDIALS_IDA_LIB NAMES sundials_ida PATHS ${IDA_DIR} NO_DEFAULT_PATH)
#find_library(SUNDIALS_SUNLINSOLLAPACKBAND_LIB NAMES sundials_sunlinsollapackband PATHS ${IDA_DIR} NO_DEFAULT_PATH)
#find_library(COOLPROP_LIB NAMES CoolProp PATHS ${COOLPROP_DIR}/64bit NO_DEFAULT_PATH)

# Search for source files
file(GLOB SRC ${PROJECT_SOURCE_DIR}/src/*.cpp)

# Add executable
add_executable(TORCCVODE ${SRC})

# # put the whole target into the bin folder
# add_custom_command(TARGET TORCCVODE
#         POST_BUILD
#         COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:TwoPhaseDropletCombustion> ${CMAKE_INSTALL_PREFIX}/bin/)


target_link_directories(TORCCVODE PRIVATE ${IDA_DIR} ${COOLPROP_DIR}/64bit ${FMT_DIR})

# Link required libraries
target_link_libraries(TORCCVODE
        PRIVATE
        PkgConfig::CANTERA
        PkgConfig::GSL
        Eigen3::Eigen
        sundials_nvecopenmp
        sundials_ida
        sundials_sunlinsollapackband
        sundials_sunlinsolspgmr
        sundials_cvode
        CoolProp
        fmt
        OpenMP::OpenMP_CXX)

# Include directories
target_include_directories(TORCCVODE
        PRIVATE ${PROJECT_SOURCE_DIR}/include ${IDA_INCLUDES} ${COOLPROP_INCLUDES} ${FMT_INCLUDES})

# Set the output path for installation
install(TARGETS TORCCVODE
        DESTINATION bin)
