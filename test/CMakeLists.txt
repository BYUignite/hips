set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib") # to include hips (built in this project)

############################################################################

add_executable(test.x)
target_sources(test.x PRIVATE test.cc)

#--------------------

set_property(TARGET test.x PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)  # adds LC_RPATH to executable

target_include_directories(test.x PRIVATE "${PROJECT_SOURCE_DIR}/src")
target_include_directories(test.x PRIVATE "$ENV{CONDA_PREFIX}/include")
target_link_directories(   test.x PRIVATE "$ENV{CONDA_PREFIX}/lib")

#--------------------

find_package(Catch2 3 REQUIRED)
target_link_libraries(test.x PRIVATE Catch2::Catch2WithMain hips)

#--------------------

install(TARGETS test.x DESTINATION "${PROJECT_SOURCE_DIR}/run")

############################################################################

if(REACTIONS_ENABLED)

    add_executable(test_rxn.x)
    target_sources(test_rxn.x PRIVATE test_rxn.cc)

    #--------------------

    set_property(TARGET test_rxn.x PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)  # adds LC_RPATH to executable

    #--------------------

    find_library(CANTERA NAMES cantera_shared cantera)
    find_package(Catch2 3 REQUIRED)

    target_include_directories(test_rxn.x PRIVATE "${PROJECT_SOURCE_DIR}/src")
    target_include_directories(test_rxn.x PRIVATE "$ENV{CONDA_PREFIX}/include")
    target_link_directories(   test_rxn.x PRIVATE "$ENV{CONDA_PREFIX}/lib")
    target_link_libraries(test_rxn.x PRIVATE hips Catch2::Catch2WithMain ${CANTERA})

    #--------------------

    install(TARGETS test_rxn.x DESTINATION "${PROJECT_SOURCE_DIR}/run")

endif()

############################################################################

add_executable(tdd.x)
target_sources(tdd.x PRIVATE tdd.cc)

#--------------------

set_property(TARGET tdd.x PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)  # adds LC_RPATH to executable

target_include_directories(tdd.x PRIVATE "${PROJECT_SOURCE_DIR}/src")
target_include_directories(tdd.x PRIVATE "$ENV{CONDA_PREFIX}/include")
target_link_directories(   tdd.x PRIVATE "$ENV{CONDA_PREFIX}/lib")

#--------------------

find_package(Catch2 3 REQUIRED)
target_link_libraries(tdd.x PRIVATE Catch2::Catch2WithMain hips)

#--------------------

install(TARGETS tdd.x DESTINATION "${PROJECT_SOURCE_DIR}/run")

############################################################################
