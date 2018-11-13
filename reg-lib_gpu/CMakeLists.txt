
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

IF(USE_CUDA)
        BUILD_INSTALL_LIBRARY(_reg_cudaCommon ${LIB_TYPE})
        BUILD_INSTALL_LIBRARY(_reg_resampling_gpu ${LIB_TYPE})
        BUILD_INSTALL_LIBRARY(_reg_tools_gpu ${LIB_TYPE})
        BUILD_INSTALL_LIBRARY(_reg_affineTransformation_gpu ${LIB_TYPE})
        BUILD_INSTALL_LIBRARY(_reg_mutualinformation_gpu ${LIB_TYPE})
        BUILD_INSTALL_LIBRARY(_reg_bspline_gpu ${LIB_TYPE})
        BUILD_INSTALL_LIBRARY(_reg_blockMatching_gpu ${LIB_TYPE})
ENDIF(USE_CUDA)

#-----------------------------------------------------------------------------
INSTALL(FILES _reg_blocksize_gpu.h DESTINATION include)

#-----------------------------------------------------------------------------