
FIND_PACKAGE(Doxygen)
if (NOT DOXYGEN_FOUND)
    message(FATAL_ERROR "Doxygen is needed to build the documentation. Please install it correctly")
ENDIF()
configure_file(Doxyfile.in ${PROJECT_BINARY_DIR}/documentation/Doxygen/Doxyfile  @ONLY IMMEDIATE)
add_custom_target(doxygen ALL COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/documentation/Doxygen/Doxyfile SOURCES ${PROJECT_BINARY_DIR}/documentation/Doxygen/Doxyfile)
INSTALL(DIRECTORY ${PROJECT_BINARY_DIR}/documentation/Doxygen/html DESTINATION ${DOC_INSTALL_DIR}/Doxygen)
