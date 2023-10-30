##################### HiPS LIBRARY #####################

add_library(hips)

target_sources(hips PRIVATE hips.cc                 hips.h 
                            batchReactor_cvode.cc   batchReactor_cvode.h   integrator_cvode.h
                            batchReactor_cantera.cc batchReactor_cantera.h
                            batchReactor.h
              )

# option(USE_CUSTOM_INTEGRATOR "Use a custom integrator" OFF){{{
# set(CUSTOM_INTEGRATOR_PATH "" CACHE STRING "Path to custom integrator library")

# if(USE_CUSTOM_INTEGRATOR)
#     target_sources(hips PRIVATE ${CUSTOM_INTEGRATOR_PATH}/customIntegrator.cc)
#     target_include_directories(hips PRIVATE ${CUSTOM_INTEGRATOR_PATH})
# else()
#     target_sources(hips PRIVATE cvodeDriver.cc cvodeDriver.h)
# endif()}}}

#----------

set_target_properties(hips PROPERTIES PUBLIC_HEADER "hips.h")

if(APPLE)
    find_library(ACCELERATE Accelerate)
    find_library(YAML yaml-cpp)
endif()

find_library(CVODES sundials_cvodes)
find_library(CANTERA cantera)
find_package(Python COMPONENTS Interpreter Development)

#----------

target_include_directories(hips PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(hips PRIVATE "/usr/local/include")
target_link_directories(   hips PRIVATE "/usr/local/lib")

target_link_libraries(hips ${CANTERA})
target_link_libraries(hips ${CANTERA} ${YAML} ${CVODES} ${Python_LIBRARIES} ${Python} ${ACCELERATE})

############# INSTALLATION ################

include(GNUInstallDirs)

install(TARGETS hips
        RUNTIME DESTINATION       "${CMAKE_INSTALL_BINDIR}"
        LIBRARY DESTINATION       "${CMAKE_INSTALL_LIBDIR}"
        ARCHIVE DESTINATION       "${CMAKE_INSTALL_LIBDIR}"
        PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}" 
       )

