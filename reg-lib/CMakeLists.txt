
#-----------------------------------------------------------------------------
SET(NAME _reg_affineTransformation)
ADD_LIBRARY(${NAME} ${LIB_TYPE} ${NAME}.h ${NAME}.cpp)
TARGET_LINK_LIBRARIES(${NAME} reg_nifti)
INSTALL(TARGETS ${NAME} DESTINATION ${LIB_INSTALL_DIR})
INSTALL(FILES ${NAME}.h DESTINATION ${INC_INSTALL_DIR})
#-----------------------------------------------------------------------------
SET(NAME _reg_resampling)
ADD_LIBRARY(${NAME} ${LIB_TYPE} ${NAME}.h ${NAME}.cpp)
TARGET_LINK_LIBRARIES(${NAME} reg_nifti _reg_affineTransformation)
INSTALL(TARGETS ${NAME} DESTINATION ${LIB_INSTALL_DIR})
INSTALL(FILES ${NAME}.h DESTINATION ${INC_INSTALL_DIR})
#-----------------------------------------------------------------------------
SET(NAME _reg_tools)
ADD_LIBRARY(${NAME} ${LIB_TYPE} ${NAME}.h ${NAME}.cpp)
TARGET_LINK_LIBRARIES(${NAME} reg_nifti)
INSTALL(TARGETS ${NAME} DESTINATION ${LIB_INSTALL_DIR})
INSTALL(FILES ${NAME}.h DESTINATION ${INC_INSTALL_DIR})
#-----------------------------------------------------------------------------
SET(NAME _reg_bspline)
ADD_LIBRARY(${NAME} ${LIB_TYPE} ${NAME}.h ${NAME}.cpp)
TARGET_LINK_LIBRARIES(${NAME} reg_nifti _reg_affineTransformation)
INSTALL(TARGETS ${NAME} DESTINATION ${LIB_INSTALL_DIR})
INSTALL(FILES ${NAME}.h DESTINATION ${INC_INSTALL_DIR})
#-----------------------------------------------------------------------------
SET(NAME _reg_mutualinformation)
ADD_LIBRARY(${NAME} ${LIB_TYPE} ${NAME}.h ${NAME}.cpp)
TARGET_LINK_LIBRARIES(${NAME} reg_nifti)
INSTALL(TARGETS ${NAME} DESTINATION ${LIB_INSTALL_DIR})
INSTALL(FILES ${NAME}.h DESTINATION ${INC_INSTALL_DIR})
#-----------------------------------------------------------------------------
SET(NAME _reg_ssd)
ADD_LIBRARY(${NAME} ${LIB_TYPE} ${NAME}.h ${NAME}.cpp)
TARGET_LINK_LIBRARIES(${NAME} reg_nifti)
INSTALL(TARGETS ${NAME} DESTINATION ${LIB_INSTALL_DIR})
INSTALL(FILES ${NAME}.h DESTINATION ${INC_INSTALL_DIR})
#-----------------------------------------------------------------------------
SET(NAME _reg_blockMatching)
ADD_LIBRARY(${NAME} ${LIB_TYPE} ${NAME}.h ${NAME}.cpp)
TARGET_LINK_LIBRARIES(${NAME} reg_nifti)
INSTALL(TARGETS ${NAME} DESTINATION ${LIB_INSTALL_DIR})
INSTALL(FILES ${NAME}.h DESTINATION ${INC_INSTALL_DIR})
#-----------------------------------------------------------------------------

MACRO(BUILD_INSTALL_LIBRARY LIBRARY_NAME LIBRARY_TYPE LIBRARY_USE_GPU)
    ADD_LIBRARY(${LIBRARY_NAME} ${LIBRARY_TYPE} ${LIBRARY_NAME}.h ${LIBRARY_NAME}.cpp)
    IF(${LIBRARY_USE_GPU} MATCHES ON)
        TARGET_LINK_LIBRARIES(${LIBRARY_NAME} ${NIFTY_REG_LIBRARIES} ${NIFTY_REG_LIBRARIES_GPU} ${CUDA_LIBRARIES} ${CUDA_cufft_LIBRARY} _niftyrec_memory _et_common)
        CUDA_ADD_CUFFT_TO_TARGET(${LIBRARY_NAME})
        CUDA_ADD_CUBLAS_TO_TARGET(${LIBRARY_NAME})
    ELSE()
        TARGET_LINK_LIBRARIES(${LIBRARY_NAME} ${NIFTY_REG_LIBRARIES} _niftyrec_memory _et_common)
    ENDIF() 
    INSTALL(TARGETS ${LIBRARY_NAME} DESTINATION ${LIB_INSTALL_DIR})
    INSTALL(FILES ${LIBRARY_NAME}.h DESTINATION ${INC_INSTALL_DIR})
ENDMACRO(BUILD_INSTALL_LIBRARY)

#-----------------------------------------------------------------------------

SET(LIBRARY_NAME _reg)
BUILD_INSTALL_LIBRARY(${LIBRARY_NAME} ${LIB_TYPE} ${NiftyRec_USE_CUDA})
SET(LIBRARY_NAME _reg_array_interface)
BUILD_INSTALL_LIBRARY(${LIBRARY_NAME} ${LIB_TYPE} ${NiftyRec_USE_CUDA})
TARGET_LINK_LIBRARIES(${LIBRARY_NAME} _reg)


