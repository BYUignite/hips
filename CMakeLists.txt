
cmake_minimum_required(VERSION 3.15)

project(hips VERSION 0.1.0 LANGUAGES CXX)

################### PROJECT SETUP ##########################

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)       # Release or Debug
endif()

# ---- Set C++ language standard ----

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS NO)


set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}/installed")

# ---- Build options ----

option(REACTIONS_ENABLED     "Enable reactions"           OFF)
option(HIPS_BUILD_EXAMPLES   "Build hips examples"        ON)
option(HIPS_BUILD_DOCS       "Build hips documentation"   OFF)
option(HIPS_BUILD_TESTS      "Build hips tests"           ON)

if(REACTIONS_ENABLED)
    add_compile_definitions(REACTIONS_ENABLED)
endif()

################### MAIN PROJECT TARGETS #####################

add_subdirectory(src)

################### OPTIONAL PROJECT TARGETS #################

# -------- EXAMPLES --------

if (HIPS_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

# -------- DOCUMENTATION ----------

if (HIPS_BUILD_DOCS )
    add_subdirectory(docs)
endif()

# -------- TESTS ----------

if (HIPS_BUILD_TESTS)
    add_subdirectory(test)
endif()
