
#  CMakeLists.txt
#  
#  NiftyRec
#  Stefano Pedemonte, May 2013.
#  CMIC - Centre for Medical Image Computing 
#  UCL - University College London. 
#  Released under BSD licence, see LICENSE.txt 

#-----------------------------------------------------------------------------
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/lib)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/teem/include/)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/niftyrec_version.in.h ${CMAKE_CURRENT_SOURCE_DIR}/niftyrec_version.h)

IF(NiftyRec_SUPPORT_NRRD)
    SET(LIB_TEEM teem )
ENDIF(NiftyRec_SUPPORT_NRRD)
#CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/rec_spect.in.xml ${CMAKE_CURRENT_SOURCE_DIR}/rec_spect.xml)
IF(NiftyRec_USE_CUDA)
    CUDA_ADD_EXECUTABLE(rec_spect rec_spect.cpp niftyrec_version.h )
    TARGET_LINK_LIBRARIES(rec_spect _reg_resampling _reg_affineTransformation _reg_tools ${LIB_TEEM} ${ZLIB} reg_nifti _reg_resampling_gpu  _reg_affineTransformation_gpu _reg_tools_gpu _reg_cudaCommon _et _et_array_interface )
ELSE(NiftyRec_USE_CUDA)
    ADD_EXECUTABLE(rec_spect rec_spect.cpp rec_spect.h niftyrec_version.h )
    TARGET_LINK_LIBRARIES(rec_spect _reg_resampling _reg_affineTransformation _reg_tools reg_nifti _et _et_array_interface ${LIB_TEEM} ${ZLIB})
ENDIF(NiftyRec_USE_CUDA)
INSTALL(TARGETS rec_spect DESTINATION ${BIN_INSTALL_DIR}/command_line_tools)
#INSTALL(FILES rec_spect.xml DESTINATION ${BIN_INSTALL_DIR})

#-----------------------------------------------------------------------------
