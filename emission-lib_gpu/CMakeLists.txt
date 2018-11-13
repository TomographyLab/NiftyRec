
#  CMakeLists.txt
#  
#  NiftyRec
#  Stefano Pedemonte, May 2012.
#  CMIC - Centre for Medical Image Computing 
#  UCL - University College London. 
#  Released under BSD licence, see LICENSE.txt 

#-----------------------------------------------------------------------------
MACRO(BUILD_INSTALL_LIBRARY LIBRARY_NAME LIBRARY_TYPE)
#        IF(${LIBRARY_TYPE} MATCHES SHARED)
#                CUDA_ADD_LIBRARY(${LIBRARY_NAME} SHARED ${LIBRARY_NAME}.h ${LIBRARY_NAME}.cu ${LIBRARY_NAME}_kernels.cu)
#                INSTALL(TARGETS ${LIBRARY_NAME} LIBRARY DESTINATION ${LIB_INSTALL_DIR})
#        ELSE()
                CUDA_ADD_LIBRARY(${LIBRARY_NAME} STATIC ${LIBRARY_NAME}.h ${LIBRARY_NAME}.cu ${LIBRARY_NAME}_kernels.cu)
                INSTALL(TARGETS ${LIBRARY_NAME} ARCHIVE DESTINATION ${LIB_INSTALL_DIR})
#        ENDIF()
        INSTALL(FILES ${LIBRARY_NAME}.h DESTINATION ${INC_INSTALL_DIR})
ENDMACRO(BUILD_INSTALL_LIBRARY)

#-----------------------------------------------------------------------------
IF(NiftyRec_USE_CUDA)
        FOREACH(LIBRARY_NAME ${NIFTY_REC_ET_LIBRARIES_GPU})
                BUILD_INSTALL_LIBRARY(${LIBRARY_NAME} ${LIB_TYPE})
        ENDFOREACH(LIBRARY_NAME)
ENDIF(NiftyRec_USE_CUDA)

#-----------------------------------------------------------------------------


