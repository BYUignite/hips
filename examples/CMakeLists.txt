
####################### EXAMPLES #######################

if(APPLE)
    find_library(ACCELERATE Accelerate)
    find_library(YAML yaml-cpp)
endif()

find_library(CVODES sundials_cvodes)
find_library(CANTERA cantera)
find_package(Python COMPONENTS Interpreter Development)

#----------

add_executable(ex_1.x ex_1.cc)
target_include_directories(ex_1.x PRIVATE "${CMAKE_SOURCE_DIR}/src" "/usr/local/include")
target_link_libraries(ex_1.x PRIVATE hips ${CANTERA} ${YAML} ${CVODES} ${Python} ${ACCELERATE})
install(TARGETS ex_1.x DESTINATION "${CMAKE_SOURCE_DIR}/run")

#----------

#add_executable(ex_2.x ex_2.cc)
#target_include_directories(ex_2.x PRIVATE "${CMAKE_SOURCE_DIR}/src" "/usr/local/include")
#target_link_libraries(ex_2.x PRIVATE hips ${CANTERA} ${YAML} ${CVODES} ${Python} ${ACCELERATE})
#install(TARGETS ex_2.x DESTINATION "${CMAKE_SOURCE_DIR}/run")

##----------

#add_executable(ex_3.x ex_3.cc)
#target_include_directories(ex_3.x PRIVATE "${CMAKE_SOURCE_DIR}/src" "/usr/local/include")
#target_link_libraries(ex_3.x PRIVATE hips ${CANTERA} ${YAML} ${CVODES} ${Python} ${ACCELERATE})
#install(TARGETS ex_3.x DESTINATION "${CMAKE_SOURCE_DIR}/run")

