#-----------------------------------------------------------------------------
IF(USE_ZLIB)
    ADD_DEFINITIONS(-DHAVE_ZLIB)
ENDIF(USE_ZLIB)
SET(NAME reg_nifti)
ADD_LIBRARY(${NAME} ${LIB_TYPE} nifti1.h nifti1_io.c nifti1_io.h znzlib.c znzlib.h)
INSTALL(TARGETS ${NAME} DESTINATION ${LIB_INSTALL_DIR})
INSTALL(FILES nifti1.h nifti1_io.h znzlib.h DESTINATION ${INC_INSTALL_DIR})
#-----------------------------------------------------------------------------
