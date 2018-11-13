
#  CMakeLists.txt
#  
#  NiftyRec
#  Stefano Pedemonte, Oct. 2012.
#  CMIC - Centre for Medical Image Computing 
#  UCL - University College London. 
#  Released under BSD licence, see LICENSE.txt 

#-----------------------------------------------------------------------------

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

MACRO(BUILD_INSTALL_SEG_MEX MEX_NAME)
    ADD_LIBRARY(${MEX_NAME} MODULE _${MEX_NAME}.cpp)
    IF(USE_GPU)
        TARGET_LINK_LIBRARIES(${MEX_NAME} ${MATLAB_MEX_LIBRARY} ${MATLAB_MX_LIBRARY} ${NIFTY_REG_LIBRARIES} _seg_EM _seg_tools _seg_MRF _seg_FMM _seg_BiasCorrection _seg_array_interface ${ZLIB})
        CUDA_ADD_CUFFT_TO_TARGET(${MEX_NAME})
        CUDA_ADD_CUBLAS_TO_TARGET(${MEX_NAME})
    ELSE(USE_GPU)
        TARGET_LINK_LIBRARIES(${MEX_NAME} ${MATLAB_MEX_LIBRARY} ${MATLAB_MX_LIBRARY} ${NIFTY_REG_LIBRARIES} _seg_EM _seg_tools _seg_MRF _seg_FMM _seg_BiasCorrection _seg_array_interface ${ZLIB})
    ENDIF(USE_GPU)
    MEXIFY(${MEX_NAME})
    INSTALL(TARGETS ${MEX_NAME} LIBRARY DESTINATION ${MEX_INSTALL_DIR_NIFTYSEG})
ENDMACRO(BUILD_INSTALL_SEG_MEX)


#-----------------------------------------------------------------------------

IF(BUILD_MATLAB)
        FOREACH(MEX_NAME ${NIFTY_REC_ET_MEX})
                BUILD_INSTALL_ET_MEX(${MEX_NAME} ${USE_CUDA})
        ENDFOREACH(MEX_NAME)
ENDIF(BUILD_MATLAB)

#-----------------------------------------------------------------------------

BUILD_INSTALL_SEG_MEX(seg_initialise)
BUILD_INSTALL_SEG_MEX(seg_set_MRF_strength)
BUILD_INSTALL_SEG_MEX(seg_get_MRF_strength)
BUILD_INSTALL_SEG_MEX(seg_set_regularisation_covariance)
BUILD_INSTALL_SEG_MEX(seg_get_regularisation_covariance)
BUILD_INSTALL_SEG_MEX(seg_step)
BUILD_INSTALL_SEG_MEX(seg_step_Gaussian)
BUILD_INSTALL_SEG_MEX(seg_step_Expectation)
BUILD_INSTALL_SEG_MEX(seg_step_Maximization)
BUILD_INSTALL_SEG_MEX(seg_step_MRF)
BUILD_INSTALL_SEG_MEX(seg_step_BiasField)
BUILD_INSTALL_SEG_MEX(seg_step_PriorWeight)
BUILD_INSTALL_SEG_MEX(seg_save_segmentation)
BUILD_INSTALL_SEG_MEX(seg_get_segmentation)
BUILD_INSTALL_SEG_MEX(seg_get_mask)
BUILD_INSTALL_SEG_MEX(seg_save_mask)
BUILD_INSTALL_SEG_MEX(seg_set_mask)
BUILD_INSTALL_SEG_MEX(seg_set_priors)
BUILD_INSTALL_SEG_MEX(seg_terminate)
BUILD_INSTALL_SEG_MEX(seg_set_biasfield_parameters)
BUILD_INSTALL_SEG_MEX(seg_get_biasfield_parameters)
BUILD_INSTALL_SEG_MEX(seg_get_biasfield)
BUILD_INSTALL_SEG_MEX(seg_save_biasfield)
BUILD_INSTALL_SEG_MEX(seg_get_image)
BUILD_INSTALL_SEG_MEX(seg_get_priors)
BUILD_INSTALL_SEG_MEX(seg_save_priors)
BUILD_INSTALL_SEG_MEX(seg_set_input_image)
BUILD_INSTALL_SEG_MEX(seg_set_segmentation)
BUILD_INSTALL_SEG_MEX(seg_get_loglikelihood)
BUILD_INSTALL_SEG_MEX(seg_get_mean_variance)
BUILD_INSTALL_SEG_MEX(seg_set_mean)
BUILD_INSTALL_SEG_MEX(seg_set_variance)
BUILD_INSTALL_SEG_MEX(seg_erode)

#-----------------------------------------------------------------------------

