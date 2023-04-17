
add_executable(hips.x)
target_sources(hips.x PRIVATE driver.cc cvodeDriver.cc cvodeDriver.h hips.cc hips.h)

find_library(CVODES sundials_cvodes)
find_library(CANTERA cantera)
find_package(Python COMPONENTS Interpreter Development)
if(APPLE)
    find_library(ACCELERATE Accelerate)
endif()

target_include_directories(hips.x PRIVATE "/usr/local/include")
target_link_libraries(hips.x ${CANTERA} ${CVODES} ${Python_LIBRARIES} ${ACCELERATE})

install(TARGETS hips.x DESTINATION "${PROJECT_SOURCE_DIR}/run")
