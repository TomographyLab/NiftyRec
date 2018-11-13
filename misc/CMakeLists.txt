
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/count_lines.py ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
EXECUTE_PROCESS(COMMAND python ${CMAKE_CURRENT_BINARY_DIR}/count_lines.py ${CMAKE_SOURCE_DIR} /build /teem /nifti /reg-lib /reg-lib_gpu /seg-lib /tclap)

