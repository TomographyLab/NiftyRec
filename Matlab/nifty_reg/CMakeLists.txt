
SUBDIRS(doc)
SUBDIRS(m)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/lib)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/mex)

INCLUDE_DIRECTORIES(${MATLAB_INCLUDE_DIRS})

MACRO(MEXIFY TARGET)
    SET_TARGET_PROPERTIES(${TARGET} PROPERTIES PREFIX "" LINKER_LANGUAGE CXX)
    IF(WIN32)
        IF(CMAKE_CL_64)
            #MESSAGE("Win 64")
            SET_TARGET_PROPERTIES(${TARGET} PROPERTIES SUFFIX .mexw64 LINK_FLAGS /export:mexFunction)
        ELSE()
            #MESSAGE("Win 32")
            SET_TARGET_PROPERTIES(${TARGET} PROPERTIES SUFFIX .mexw32 LINK_FLAGS /export:mexFunction)
        ENDIF()
    ELSEIF(APPLE)
        IF(CMAKE_SIZEOF_VOID_P MATCHES "8")
            #MESSAGE("Mac 64")
            SET_TARGET_PROPERTIES(${TARGET} PROPERTIES SUFFIX .mexmaci64)
        ELSE()
            #MESSAGE("Mac 32")
            SET_TARGET_PROPERTIES(${TARGET} PROPERTIES SUFFIX .mexmaci)
        ENDIF()
    ELSEIF(UNIX)
        IF(CMAKE_SIZEOF_VOID_P MATCHES "8")
            #MESSAGE("Linux 64")
            SET_TARGET_PROPERTIES(${TARGET} PROPERTIES SUFFIX .mexa64)
        ELSE()
            #MESSAGE("Linux 32")
            SET_TARGET_PROPERTIES(${TARGET} PROPERTIES SUFFIX .mexglx)
        ENDIF ()
    ENDIF()
ENDMACRO(MEXIFY)

MACRO(BUILD_INSTALL_REG_MEX MEX_NAME)
    ADD_LIBRARY(${MEX_NAME} MODULE _${MEX_NAME}.cpp)
    IF(USE_GPU)
        TARGET_LINK_LIBRARIES(${MEX_NAME} ${MATLAB_MEX_LIBRARY} ${MATLAB_MX_LIBRARY} _reg _reg_array_interface ${NIFTY_REG_LIBRARIES} ${NIFTY_REC_REG_INTERF_LIBRARY} ${ZLIB} ${CUDA_LIBRARIES})
        CUDA_ADD_CUFFT_TO_TARGET(${MEX_NAME})
        CUDA_ADD_CUBLAS_TO_TARGET(${MEX_NAME})
    ELSE(USE_GPU)
        TARGET_LINK_LIBRARIES(${MEX_NAME} ${MATLAB_MEX_LIBRARY} ${MATLAB_MX_LIBRARY} _reg _reg_array_interface ${NIFTY_REG_LIBRARIES} ${NIFTY_REC_REG_INTERF_LIBRARY} ${ZLIB})
    ENDIF(USE_GPU)
    MEXIFY(${MEX_NAME})
    INSTALL(TARGETS ${MEX_NAME} LIBRARY DESTINATION ${MEX_INSTALL_DIR_NIFTYREG})
ENDMACRO(BUILD_INSTALL_REG_MEX)

#-----------------------------------------------------------------------------

BUILD_INSTALL_REG_MEX(reg_control_points_from_affine_mex)
BUILD_INSTALL_REG_MEX(reg_gradient_NMI_nodes_mex)
BUILD_INSTALL_REG_MEX(reg_gradient_voxel_to_nodes_mex)
BUILD_INSTALL_REG_MEX(reg_resample_spline_mex)
BUILD_INSTALL_REG_MEX(reg_image_gradient_mex)
BUILD_INSTALL_REG_MEX(reg_ssd_gradient_mex)
BUILD_INSTALL_REG_MEX(reg_gaussian_smooth_mex)
BUILD_INSTALL_REG_MEX(reg_scale_amplitude_mex)
BUILD_INSTALL_REG_MEX(reg_gradient_bending_energy_mex)
BUILD_INSTALL_REG_MEX(reg_gradient_jacobian_determinant_mex)


#-----------------------------------------------------------------------------



