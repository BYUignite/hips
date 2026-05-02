find_program(Python3_EXECUTABLE NAMES python python3)
#find_package(Python3 COMPONENTS Interpreter)

if (Python3_EXECUTABLE AND BUILD_PYTHON_INTERFACE)

    set(PYHIPS_PACKAGE_NAME "pyhips")
    set(PYHIPS_PACKAGE_VERSION "${PROJECT_VERSION}")

    set(PYHIPS_EXTENSION_SOURCES
        "'${CMAKE_CURRENT_SOURCE_DIR}/pyhips.pyx'"
        "'${CMAKE_CURRENT_SOURCE_DIR}/../c++/hips.cc'"
    )

    if(REACTIONS_ENABLED)
        list(APPEND PYHIPS_EXTENSION_SOURCES
            "'${CMAKE_CURRENT_SOURCE_DIR}/../c++/batchReactor_cvode.cc'"
            "'${CMAKE_CURRENT_SOURCE_DIR}/../c++/batchReactor_cantera.cc'"
        )
    endif()

    get_target_property(HIPS_INCLUDE_DIRS hips INCLUDE_DIRECTORIES)
    if(NOT HIPS_INCLUDE_DIRS)
        set(HIPS_INCLUDE_DIRS)
    endif()

    get_target_property(HIPS_LINK_DIRECTORIES hips LINK_DIRECTORIES)
    if(NOT HIPS_LINK_DIRECTORIES)
        set(HIPS_LINK_DIRECTORIES)
    endif()

    get_target_property(HIPS_LINK_LIBRARIES hips LINK_LIBRARIES)
    if(NOT HIPS_LINK_LIBRARIES)
        set(HIPS_LINK_LIBRARIES)
    endif()

    get_target_property(HIPS_COMPILE_DEFINITIONS hips COMPILE_DEFINITIONS)
    if(NOT HIPS_COMPILE_DEFINITIONS)
        set(HIPS_COMPILE_DEFINITIONS)
    endif()

    if(REACTIONS_ENABLED)
        list(APPEND HIPS_COMPILE_DEFINITIONS REACTIONS_ENABLED)
    endif()
    list(REMOVE_DUPLICATES HIPS_COMPILE_DEFINITIONS)

    list(JOIN PYHIPS_EXTENSION_SOURCES ",\n                " PYHIPS_EXTENSION_SOURCES_PY)

    set(PYCGAS_EXTENSION_BLOCK "")
    if(REACTIONS_ENABLED)
        set(PYCGAS_EXTENSION_SOURCES
            "'${CMAKE_CURRENT_SOURCE_DIR}/pycgas.pyx'"
        )
        list(JOIN PYCGAS_EXTENSION_SOURCES ",\n                " PYCGAS_EXTENSION_SOURCES_PY)
        set(PYCGAS_EXTENSION_BLOCK "
    Extension('pycgas',
              sources=[${PYCGAS_EXTENSION_SOURCES_PY}],
              language='c++',
              include_dirs=include_dirs,
              library_dirs=library_dirs,
              define_macros=define_macros,
              extra_link_args=link_args,
              extra_compile_args=['-std=c++20', '-O3']),")
    endif()

    set(PYHIPS_INCLUDE_DIRS "    '${CMAKE_CURRENT_SOURCE_DIR}/../c++',\n")
    if(REACTIONS_ENABLED)
        string(APPEND PYHIPS_INCLUDE_DIRS "    '${CMAKE_CURRENT_SOURCE_DIR}/../c++/test',\n")
    endif()
    foreach(dir IN LISTS HIPS_INCLUDE_DIRS)
        if(dir MATCHES "\\$<")
            continue()
        endif()
        string(APPEND PYHIPS_INCLUDE_DIRS "    '${dir}',\n")
    endforeach()

    set(PYHIPS_LIBRARY_DIRS "")
    foreach(dir IN LISTS HIPS_LINK_DIRECTORIES)
        string(APPEND PYHIPS_LIBRARY_DIRS "    '${dir}',\n")
    endforeach()

    set(PYHIPS_LINK_ARGS "")
    foreach(lib IN LISTS HIPS_LINK_LIBRARIES)
        string(APPEND PYHIPS_LINK_ARGS "    '${lib}',\n")
    endforeach()

    set(PYHIPS_DEFINE_MACROS "")
    list(REMOVE_DUPLICATES HIPS_COMPILE_DEFINITIONS)
    foreach(def IN LISTS HIPS_COMPILE_DEFINITIONS)
        string(APPEND PYHIPS_DEFINE_MACROS "    ('${def}', None),\n")
    endforeach()

    set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
    set(SETUP_PY    "${CMAKE_CURRENT_SOURCE_DIR}/setup.py")
    configure_file(${SETUP_PY_IN} ${SETUP_PY})

    set(OUTPUT      "${CMAKE_BINARY_DIR}/python")

    add_custom_command(OUTPUT ${OUTPUT}
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                       COMMAND ${Python3_EXECUTABLE} ${SETUP_PY} build_ext --inplace
                       COMMAND cmake -E touch ${OUTPUT})

    add_custom_target(python ALL DEPENDS ${OUTPUT})

    install(CODE "execute_process(COMMAND ${Python3_EXECUTABLE} ${SETUP_PY} install --prefix=${CMAKE_INSTALL_PREFIX})")

elseif(NOT Python3_EXECUTABLE AND BUILD_PYTHON_INTERFACE)

    message(FATAL_ERROR "Python not found")

endif()
