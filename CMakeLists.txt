cmake_minimum_required(VERSION 3.19)
project("invertible-cdf" VERSION 1.0
  DESCRIPTION "A header only cmake/c++ library exposing <InvertibleCDF> implementation"
  HOMEPAGE_URL "https://github.com/DominikHorn/invertible-cdf"
  LANGUAGES CXX)

# Declare library & directories to include. See
# http://mariobadr.com/creating-a-header-only-library-with-cmake.html for more
# info/install instructions
add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(${PROJECT_NAME} INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

# ==== Dependencies ====
include(${PROJECT_SOURCE_DIR}/thirdparty/learned_hashing.cmake)
target_link_libraries(${PROJECT_NAME} INTERFACE ${LEARNED_HASHING_LIBRARY})

# Make IDE friendly
target_sources(${PROJECT_NAME} INTERFACE invertible-cdf.hpp include/)

# Benchmark and test code. Only add_subdirectory() if
#   this code is used standalone, i.e., not as a FetchContent()
#   or other dependency. Standalone detection uses `hasParent`
get_directory_property(hasParent PARENT_DIRECTORY)
if (NOT hasParent)
  add_subdirectory(src)
endif ()
