cmake_minimum_required(VERSION 3.16)
project(SfSpp VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# --- Homebrew hints (macOS) -------------------------------------------------
if(APPLE)
  execute_process(COMMAND brew --prefix OUTPUT_VARIABLE HOMEBREW_PREFIX
                  OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  if(HOMEBREW_PREFIX)
    list(APPEND CMAKE_PREFIX_PATH "${HOMEBREW_PREFIX}")
    # AppleClang needs an explicit pointer at keg-only libomp
    if(EXISTS "${HOMEBREW_PREFIX}/opt/libomp")
      set(OpenMP_ROOT "${HOMEBREW_PREFIX}/opt/libomp" CACHE PATH "")
      set(OpenMP_CXX_FLAGS "-Xclang -fopenmp -I${HOMEBREW_PREFIX}/opt/libomp/include")
      set(OpenMP_CXX_LIB_NAMES "omp")
      set(OpenMP_omp_LIBRARY "${HOMEBREW_PREFIX}/opt/libomp/lib/libomp.dylib")
    endif()
  endif()
endif()

find_package(Ceres CONFIG REQUIRED)
find_package(PCL CONFIG REQUIRED COMPONENTS common io features visualization)
find_package(OpenMP)

file(GLOB PROJ_SRC
  "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/class/*.cpp")

add_executable(SfSpp ${PROJ_SRC})

target_include_directories(SfSpp
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/class
    ${PCL_INCLUDE_DIRS})

target_link_libraries(SfSpp
  PRIVATE
    Ceres::ceres
    ${PCL_LIBRARIES})

target_compile_definitions(SfSpp PRIVATE ${PCL_DEFINITIONS})

if(OpenMP_CXX_FOUND)
  target_link_libraries(SfSpp PRIVATE OpenMP::OpenMP_CXX)
else()
  message(WARNING "OpenMP not found - building single threaded")
endif()

if(MSVC)
  target_compile_options(SfSpp PRIVATE /bigobj)
  target_compile_definitions(SfSpp PRIVATE _CRT_SECURE_NO_WARNINGS)
else()
  target_compile_options(SfSpp PRIVATE -Wno-deprecated-declarations)
endif()

# --- extract_axis -----------------------------------------------------------
# Replaces the MATLAB-only axis stage of the preprocessing pipeline; see
# tools/extract_axis.cpp. Reuses this repository's PotSAC implementation, so it
# only needs the data structures and the axis code, not the reassembly.
add_executable(extract_axis
  ${CMAKE_CURRENT_SOURCE_DIR}/tools/extract_axis.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/class/axis_estimation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/class/data_structure.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/class/KDTree.cpp)

target_include_directories(extract_axis
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/class
    ${PCL_INCLUDE_DIRS})

target_link_libraries(extract_axis PRIVATE Ceres::ceres ${PCL_LIBRARIES})
target_compile_definitions(extract_axis PRIVATE ${PCL_DEFINITIONS})
if(OpenMP_CXX_FOUND)
  target_link_libraries(extract_axis PRIVATE OpenMP::OpenMP_CXX)
endif()
if(NOT MSVC)
  target_compile_options(extract_axis PRIVATE -Wno-deprecated-declarations)
endif()

# --- SfSpp_ctrl -------------------------------------------------------------
# Same sources, POT_CTRL selected at compile time: the authors' Pot A driven
# through our own preprocessing chain, scored against their ground truth. Built
# as a separate binary so a control run and a working run can proceed together.
add_executable(SfSpp_ctrl ${PROJ_SRC})
target_include_directories(SfSpp_ctrl PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/class ${PCL_INCLUDE_DIRS})
target_link_libraries(SfSpp_ctrl PRIVATE Ceres::ceres ${PCL_LIBRARIES})
target_compile_definitions(SfSpp_ctrl PRIVATE ${PCL_DEFINITIONS} POT_CTRL)
if(OpenMP_CXX_FOUND)
  target_link_libraries(SfSpp_ctrl PRIVATE OpenMP::OpenMP_CXX)
endif()
if(NOT MSVC)
  target_compile_options(SfSpp_ctrl PRIVATE -Wno-deprecated-declarations)
endif()
