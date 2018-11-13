
#  CMakeLists.txt
#  
#  NiftyRec
#  Stefano Pedemonte, May 2012.
#  CMIC - Centre for Medical Image Computing 
#  UCL - University College London. 
#  Released under BSD licence, see LICENSE.txt 

#-----------------------------------------------------------------------------
MACRO(BUILD_INSTALL_LIBRARY LIBRARY_NAME LIBRARY_TYPE LIBRARY_USE_GPU)
    ADD_LIBRARY(${LIBRARY_NAME} ${LIBRARY_TYPE} ${LIBRARY_NAME}.h ${LIBRARY_NAME}.cpp)
    IF(${LIBRARY_USE_GPU} MATCHES ON)
        TARGET_LINK_LIBRARIES(${LIBRARY_NAME} ${CUDA_LIBRARIES} ${CUDA_cufft_LIBRARY})
        CUDA_ADD_CUFFT_TO_TARGET(${LIBRARY_NAME})
        CUDA_ADD_CUBLAS_TO_TARGET(${LIBRARY_NAME})
    ELSE()
        #TARGET_LINK_LIBRARIES(${LIBRARY_NAME})
    ENDIF()
    INSTALL(TARGETS ${LIBRARY_NAME} DESTINATION ${LIB_INSTALL_DIR})
    INSTALL(FILES ${LIBRARY_NAME}.h DESTINATION ${INC_INSTALL_DIR})
ENDMACRO(BUILD_INSTALL_LIBRARY)

#-----------------------------------------------------------------------------

FOREACH(LIBRARY_NAME ${NIFTY_REC_TT_LIBRARIES})
    BUILD_INSTALL_LIBRARY(${LIBRARY_NAME} ${LIB_TYPE} ${NiftyRec_USE_CUDA})
ENDFOREACH(LIBRARY_NAME)

#-----------------------------------------------------------------------------

SET(LIBRARY_NAME _tt)
BUILD_INSTALL_LIBRARY(${LIBRARY_NAME} ${LIB_TYPE} ${NiftyRec_USE_CUDA})
TARGET_LINK_LIBRARIES(${LIBRARY_NAME} ${NIFTY_REC_ET_LIBRARIES} ${NIFTY_REC_TT_LIBRARIES})
IF(NiftyRec_USE_CUDA)
    TARGET_LINK_LIBRARIES(${LIBRARY_NAME} ${NIFTY_REC_TT_LIBRARIES_GPU} ${NIFTY_REC_ET_LIBRARIES_GPU})
    CUDA_ADD_CUFFT_TO_TARGET(${LIBRARY_NAME})
    CUDA_ADD_CUBLAS_TO_TARGET(${LIBRARY_NAME})
ENDIF(NiftyRec_USE_CUDA)

SET(LIBRARY_NAME _tt_array_interface)
BUILD_INSTALL_LIBRARY(${LIBRARY_NAME} ${LIB_TYPE} ${NiftyRec_USE_CUDA})
TARGET_LINK_LIBRARIES(${LIBRARY_NAME} _tt)

#-----------------------------------------------------------------------------
