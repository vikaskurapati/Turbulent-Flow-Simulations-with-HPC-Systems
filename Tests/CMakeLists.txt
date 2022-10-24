file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS "*.cpp")

foreach(file ${SOURCES})
    get_filename_component(filename ${file} NAME_WLE)
    display_header("Creating Makefile of ${filename}")
    add_executable(${filename} ${file})
    add_test(NAME ${filename} COMMAND ${filename})
    target_link_libraries(${filename} PRIVATE ${META_PROJECT_NAME})
endforeach()

add_test(NAME Cavity2DTest
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:${META_PROJECT_NAME}-Runner> ${CMAKE_BINARY_DIR}/ExampleCases/Cavity2D.xml
)

add_test(NAME Cavity3DTest
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:${META_PROJECT_NAME}-Runner> ${CMAKE_BINARY_DIR}/ExampleCases/Cavity3D.xml
)

add_test(NAME Channel2DTest
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:${META_PROJECT_NAME}-Runner> ${CMAKE_BINARY_DIR}/ExampleCases/Channel2D.xml
)

add_test(NAME Channel3DTest
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:${META_PROJECT_NAME}-Runner> ${CMAKE_BINARY_DIR}/ExampleCases/Channel3D.xml
)

add_test(NAME ChannelBackwardFacingStep2DTest
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:${META_PROJECT_NAME}-Runner> ${CMAKE_BINARY_DIR}/ExampleCases/ChannelBackwardFacingStep2D.xml
)

add_test(NAME ChannelBackwardFacingStep3DTest
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:${META_PROJECT_NAME}-Runner> ${CMAKE_BINARY_DIR}/ExampleCases/ChannelBackwardFacingStep3D.xml
)
