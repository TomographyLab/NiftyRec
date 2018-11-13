SET(MATLAB_FOUND 0)

include(LibFindMacros)

FIND_PATH(MATLAB_INCLUDE
        NAMES "mex.h"
        PATHS "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MATLAB\\R2009b;MATLABROOT]/extern/include"
)
FIND_PATH(MATLAB_LIBS_DIR
        NAMES "libmex.lib"
        PATHS "[HKEY_LOCAL_MACHINE\\SOFTWARE\\MATLAB\\R2009b;MATLABROOT]/bin/win32"
)

SET(MATLAB_LIBRARIES "")

IF(WIN32)
    SET(LIB_SUFFIX "lib")
ELSE(WIN32)
    SET(LIB_SUFFIX "")
ENDIF(WIN32)

FIND_LIBRARY(MATLAB_MEX_LIBRARY
    NAMES "${LIB_SUFFIX}mex"
    PATHS ${MATLAB_LIBS_DIR}
    DOC "Matlab mex library"
    NO_DEFAULT_PATH
)
mark_as_advanced(MATLAB_MEX_LIBRARY)

FIND_LIBRARY(MATLAB_MX_LIBRARY
    NAMES "${LIB_SUFFIX}mx"
    PATHS ${MATLAB_LIBS_DIR}
    DOC "Matlab mx library"
    NO_DEFAULT_PATH
)
mark_as_advanced(MATLAB_MX_LIBRARY)

FIND_LIBRARY(MATLAB_ENG_LIBRARY
    NAMES "${LIB_SUFFIX}mat"
    PATHS ${MATLAB_LIBS_DIR}
    DOC "Matlab engine library"
    NO_DEFAULT_PATH
)
mark_as_advanced(MATLAB_ENG_LIBRARY)

SET(MATLAB_LIBRARY 
    "${MATLAB_MEX_LIBRARY}" 
    "${MATLAB_MX_LIBRARY}"
    "${MATLAB_ENG_LIBRARY}"
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(MATLAB_PROCESS_INCLUDES MATLAB_INCLUDE)
set(MATLAB_PROCESS_LIBS MATLAB_LIBRARY MATLAB_LIBRARIES)
libfind_process(MATLAB)
