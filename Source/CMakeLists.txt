add_library(${META_PROJECT_NAME} OBJECT)
display_header("Creating Makefile of ${META_PROJECT_NAME}")

file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "*")
list(FILTER SOURCES EXCLUDE REGEX ".*Main\\.cpp$")
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${SOURCES})
target_sources(${META_PROJECT_NAME} PRIVATE ${SOURCES})
target_include_directories(${META_PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(${META_PROJECT_NAME}
    PUBLIC
        project_options
        project_warnings
)

find_package(MPI REQUIRED)
message(STATUS "Run: ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS")
find_package(spdlog REQUIRED)
find_package(tinyxml2 REQUIRED)
target_link_system_libraries(${META_PROJECT_NAME}
    PUBLIC
        MPI::MPI_CXX
        spdlog::spdlog
        tinyxml2::tinyxml2
)

option(ENABLE_SINGLE_PRECISION "Enable single floating-point precision" OFF)
if(ENABLE_SINGLE_PRECISION)
    target_compile_definitions(${META_PROJECT_NAME} PUBLIC ENABLE_SINGLE_PRECISION)
endif()

option(ENABLE_PETSC "Enable the Portable, Extensible Toolkit for Scientific Computation (PETSc)" ON)
if(ENABLE_PETSC)
    find_package(PETSc REQUIRED)
    target_compile_definitions(${META_PROJECT_NAME} PUBLIC ENABLE_PETSC)
    target_link_system_libraries(${META_PROJECT_NAME} PUBLIC ${PETSc_LIBRARIES})
    target_include_system_directories(${META_PROJECT_NAME} PUBLIC ${PETSc_INCLUDE_DIRS})
endif()

add_executable(${META_PROJECT_NAME}-Runner Main.cpp)
target_link_libraries(${META_PROJECT_NAME}-Runner PRIVATE ${META_PROJECT_NAME})

if(MSVC)
    file(COPY ${CMAKE_SOURCE_DIR}/ExampleCases DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    configure_file(${CMAKE_SOURCE_DIR}/Tools/petsc_commandline_arg ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
else()
    file(COPY ${CMAKE_SOURCE_DIR}/ExampleCases DESTINATION ${CMAKE_BINARY_DIR})
    configure_file(${CMAKE_SOURCE_DIR}/Tools/petsc_commandline_arg ${CMAKE_BINARY_DIR} COPYONLY)
endif()
