
####################### EXAMPLES #######################

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib") # to include hips (built in this project)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)            # only needed if examples dep on an external lib

#---------- Example 1

add_executable(ex_1.x ex_1.cc)
target_include_directories(ex_1.x PRIVATE "${CMAKE_SOURCE_DIR}/src")
target_link_libraries(ex_1.x PRIVATE hips)
install(TARGETS ex_1.x DESTINATION "${PROJECT_SOURCE_DIR}/run")

#---------- Example 2

if(REACTIONS_ENABLED)

    find_library(CANTERA NAMES cantera_shared cantera)

    add_executable(ex_2.x ex_2.cc)
    target_include_directories(ex_2.x PUBLIC "$ENV{CONDA_PREFIX}/include")
    target_link_directories(ex_2.x PUBLIC "$ENV{CONDA_PREFIX}/lib")
    target_link_libraries(ex_2.x PRIVATE hips ${CANTERA})
    install(TARGETS ex_2.x DESTINATION "${PROJECT_SOURCE_DIR}/run")

    add_executable(ex_t.x ex_t.cc)
    target_include_directories(ex_t.x PUBLIC "$ENV{CONDA_PREFIX}/include")
    target_link_directories(ex_t.x PUBLIC "$ENV{CONDA_PREFIX}/lib")
    target_link_libraries(ex_t.x PRIVATE hips ${CANTERA})
    install(TARGETS ex_t.x DESTINATION "${PROJECT_SOURCE_DIR}/run")

endif()

#---------- Example 3

add_executable(ex_3.x ex_3.cc)
target_include_directories(ex_3.x PRIVATE "${CMAKE_SOURCE_DIR}/src")
target_link_libraries(ex_3.x PRIVATE hips)
install(TARGETS ex_3.x DESTINATION "${PROJECT_SOURCE_DIR}/run")

