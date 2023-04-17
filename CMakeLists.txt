################### PREAMBLE ###############################

cmake_minimum_required(VERSION 3.15)

project(hips VERSION 0.1.0 LANGUAGES CXX)

################### PROJECT SETUP ##########################

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)       # Release or Debug
endif()

# ---- set C++ language standard ----

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS NO)

# ---- installation location ----

set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}/installation")

# ---- build options ----

option(HIPS_BUILD_EXAMPLES   "Build hips examples"        ON)
option(HIPS_BUILD_DOCS       "Build hips documentation"   ON)

################### MAIN PROJECT TARGETS #####################

add_subdirectory(src)

################### OPTIONAL PROJECT TARGETS #################

# -------- EXAMPLES --------

#if (HIPS_BUILD_EXAMPLES)
#    add_subdirectory(examples)
#endif()

# -------- DOCUMENTATION --------

#if (HIPS_BUILD_DOCS)
#    add_subdirectory(docs)
#endif()
