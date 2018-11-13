
#  CMakeLists.txt
#  
#  NiftyRec
#  Stefano Pedemonte, 2009-2013. 
#  - CMIC - Centre for Medical Image Computing, UCL, London. (2009-2013)
#  - Aalto University School of Science, Helsinki. (2013-2014)
#  - A.A. Martinos Center for Biomedical Imaging, MGH, Harvard University. Boston. (2013-2014)
#  Released under BSD licence, see LICENSE.txt 


#-----------------------------------------------------------------------------
# Set default options: 

# Set default build type to Release. If -DCMAKE_BUILD_TYPE is specified on the command line, use that definition instead. 
IF(DEFINED CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of 
build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug
Release.") 
ELSE()
   SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build,
options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release.") 
ENDIF()

# Set default install prefix. If -DCMAKE_INSTALL_PREFIX is specified on the command line, use that definition instead. 
SET(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/../NiftyRec_install/ CACHE PATH "NiftyRec install prefix")

#-----------------------------------------------------------------------------
# Setup CMake project: 

PROJECT(nifty_rec)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
cmake_policy(VERSION 2.8)
INCLUDE(ExternalProject)
ENABLE_TESTING()

#-----------------------------------------------------------------------------
# NiftyRec version: 
SET(Nifty_Rec_VERSION_MAJOR 3)
SET(Nifty_Rec_VERSION_MINOR 0)
SET(Nifty_Rec_VERSION_PATCH 0)

#-----------------------------------------------------------------------------
# Build options: 
OPTION(NiftyRec_BUILD_APPS          "Build command line applications (including 3D Slicer plugins)" OFF)
OPTION(NiftyRec_BUILD_MATLAB        "Build Matlab extension modules" OFF)
IF(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    IF(NOT(${CMAKE_SYSTEM_NAME} MATCHES "Windows"))
        IF(NOT CMAKE_SYSTEM_VERSION VERSION_GREATER 11.0 OR CMAKE_SYSTEM_VERSION VERSION_EQUAL 11.0.0)
            OPTION(NiftyRec_USE_OPENMP "To use openMP for multi-CPU processing" ON)
        ENDIF()
    ENDIF()
ENDIF()
OPTION(NiftyRec_INCLUDE_DATA        "Include data sets for examples and demos" ON)
OPTION(NiftyRec_USE_CUDA            "Compile GPU implementation of the code" OFF)
OPTION(NiftyRec_SUPPORT_NRRD        "Support NRRD image file format. Enable to generate 3D Slicer plugins." OFF)
OPTION(NiftyRec_COMPILE_DOC         "Compile pdf documentation from the Latex source" OFF)
OPTION(NiftyRec_BUILD_RENDERER      "Build NiftyRec volume renderer" OFF)
OPTION(NiftyRec_VERBOSE             "Print out extra information when compiling" OFF)

IF(NiftyRec_BUILD_RENDERER)
    SET(NiftyRec_USE_CUDA ON)
    MESSAGE(STATUS "NiftyRec_USE_CUDA has been enabled in order to build the volume renderer.")
ENDIF(NiftyRec_BUILD_RENDERER)

IF(NOT(${CMAKE_SYSTEM_NAME} MATCHES "Windows"))
    SET(NiftyRec_BUILD_SHARED ON)
ELSE()
    SET(NiftyRec_BUILD_SHARED OFF)
ENDIF()
SET(NiftyRec_DO_PACKAGE ON)

IF(CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "limited configs" FORCE)
ENDIF(CMAKE_CONFIGURATION_TYPES)
SET(CMAKE_BUILD_TYPE "Release") 

#-----------------------------------------------------------------------------
# Setup OpenMP multiprocessing environment. 
FIND_PACKAGE(Threads)
IF(NiftyRec_USE_OPENMP)
    EXECUTE_PROCESS(COMMAND ${CMAKE_C_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    IF(GCC_VERSION VERSION_GREATER 4.2 OR GCC_VERSION VERSION_EQUAL 4.2)
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
    ENDIF()
    IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /openmp")
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /openmp")
    ENDIF()
ENDIF(NiftyRec_USE_OPENMP)

#-----------------------------------------------------------------------------
# Setup CUDA environment: 

SET(USE_CUDA 0)
IF(NiftyRec_USE_CUDA)
        SET(USE_CUDA 1)
        FIND_PACKAGE(CUDA REQUIRED)
        ADD_DEFINITIONS(-D_USE_CUDA)
        SET(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS};-D_USE_CUDA")
        INCLUDE_DIRECTORIES(${CUDA_TOOLKIT_INCLUDE})

        # MAC OSX 10.8 and 10.9 compatibility 
        IF(APPLE)
                EXEC_PROGRAM(uname ARGS -v  OUTPUT_VARIABLE DARWIN_VERSION)
                STRING(REGEX MATCH "[0-9]+" DARWIN_VERSION ${DARWIN_VERSION})
                MESSAGE(STATUS "DARWIN_VERSION=${DARWIN_VERSION}")
                IF(DARWIN_VERSION GREATER 12)  # Darwin version 12 is OS X 10.8
                        #SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -F/Library/Frameworks -framework CUDA")
                        #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -F/Library/Frameworks -framework CUDA")
                        SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -F/Library/Frameworks")
                        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -F/Library/Frameworks")
                ENDIF(DARWIN_VERSION GREATER 12) 
ENDIF(APPLE) 

ENDIF(NiftyRec_USE_CUDA)

#-----------------------------------------------------------------------------
# Setup the bindings for Matlab 

IF(UNIX AND NiftyRec_BUILD_MATLAB AND NOT NiftyRec_BUILD_SHARED)
        MESSAGE("Unix platform: desabling Matlab extension because libraries are not built as shared.")
        SET(NiftyRec_BUILD_MATLAB OFF)
ENDIF()

#-----------------------------------------------------------------------------
# Set build options: 

# Decide whether to build as static or shared: 
IF(NiftyRec_BUILD_SHARED)
        SET(LIB_TYPE SHARED)
        IF(UNIX)
                IF("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
                        ADD_DEFINITIONS(-fPIC)
                        IF(NOT APPLE)
                            SET(CMAKE_CC_FLAGS_RELEASE "${CMAKE_CC_FLAGS_RELEASE} -fPIC -ldl"  )
                            SET(CMAKE_CC_FLAGS_DEBUG "${CMAKE_CC_FLAGS_DEBUG} -fPIC -ldl"  )
                            SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fPIC -ldl"  )
                            SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fPIC -ldl"  )
                        ENDIF(NOT APPLE)
                ENDIF("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
        ENDIF(UNIX)
ELSE(NiftyRec_BUILD_SHARED)
        SET(LIB_TYPE STATIC)
ENDIF(NiftyRec_BUILD_SHARED)

# Print verbose information at build time: 
IF(NiftyRec_VERBOSE)
	ADD_DEFINITIONS(-D_VERBOSE)
        SET(CUDA_VERBOSE_BUILD ON)
ELSE(NiftyRec_VERBOSE)
        SET(CUDA_VERBOSE_BUILD OFF)
ENDIF(NiftyRec_VERBOSE)

# Store relative path information in each library. This avoids setting LD_LIBRARY_PATH manually: 
#SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
#SET(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/${LIB_INSTALL_DIR})

# Build Teem NRRD libraries: 
IF(NiftyRec_SUPPORT_NRRD)
        ADD_DEFINITIONS(-D_SUPPORT_NRRD)
ENDIF(NiftyRec_SUPPORT_NRRD)

# Build nifti_io libraries:
SUBDIRS(nifti)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/nifti)

# Build tclap libraries:
SUBDIRS(tclap)
IF(NiftyRec_SUPPORT_NRRD)
        SUBDIRS(teem)
ENDIF(NiftyRec_SUPPORT_NRRD)

# Build NiftyRec emission libraries:  
SET(NIFTY_REC_ET_INTERF_LIBRARY _et_array_interface)
SET(NIFTY_REC_ET_LIBRARIES _et_common _niftyrec_memory _et_line_integral _et_line_integral_attenuated _et_line_backproject_attenuated _et_line_backproject _et_accumulate _et_clear_accumulator _et_convolve2D _et_convolveSeparable2D _et_convolveFFT2D)
SET(NIFTY_REC_ET_LIBRARIES_GPU _et_line_integral_gpu _et_line_integral_attenuated_gpu _et_line_backproject_gpu _et_line_backproject_attenuated_gpu _et_clear_accumulator_gpu _et_accumulate_gpu _et_convolveSeparable2D_gpu _et_convolveFFT2D_gpu _et_attenuation_gradient_gpu _pet_line_integral_compressed_gpu _pet_line_backproject_compressed_gpu)
SET(NIFTY_REC_ET_MEX et_project_mex et_project_partial_mex et_backproject_mex et_project_backproject_mex et_convolve_mex et_rotate_mex et_irotate_mex et_list_gpus_mex et_set_gpu_mex et_affine_mex et_reorder_activity_in et_reorder_activity_out et_variance_centre_approx_mex et_fisher_grid_mex et_fisher_grid_invprojection_mex et_gradient_attenuation_mex et_isinstalled_mex et_reset_gpu_mex et_histogram_weighted_mex)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/emission-lib)
IF(NiftyRec_USE_CUDA)
    INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/emission-lib_gpu)
ENDIF(NiftyRec_USE_CUDA)

SUBDIRS(emission-lib)
IF(NiftyRec_USE_CUDA)
	SUBDIRS(emission-lib_gpu)
ENDIF(NiftyRec_USE_CUDA)

# Build NiftyRec transmission libraries: 
SET(NIFTY_REC_TT_INTERF_LIBRARY _tt_array_interface)
SET(NIFTY_REC_TT_LIBRARIES _tt_common _tt_line_project_ray_cpu _tt_line_backproject_ray_cpu)
SET(NIFTY_REC_TT_LIBRARIES_GPU _tt_line_project_ray_gpu _tt_line_backproject_ray_gpu) 
SET(NIFTY_REC_TT_MEX tt_project_ray_mex tt_backproject_ray_mex) 

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/transmission-lib)
IF(NiftyRec_USE_CUDA)
    INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/transmission-lib_gpu)
ENDIF(NiftyRec_USE_CUDA)

SUBDIRS(transmission-lib)
IF(NiftyRec_USE_CUDA)
	SUBDIRS(transmission-lib_gpu)
ENDIF(NiftyRec_USE_CUDA)

#Build NiftyReg Libraries 
SET(NIFTY_REC_REG_INTERF_LIBRARY _reg_array_interface)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/reg-lib)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/reg-lib_gpu)
SUBDIRS(reg-lib)
SUBDIRS(reg-lib_gpu)
SET(NIFTY_REG_LIBRARIES reg_nifti _reg_tools _reg_resampling _reg_affineTransformation _reg_bspline _reg_ssd _reg_mutualinformation _reg_blockMatching )
SET(NIFTY_REG_LIBRARIES_GPU _reg_cudaCommon _reg_resampling_gpu _reg_affineTransformation_gpu _reg_mutualinformation_gpu _reg_bspline_gpu _reg_tools_gpu _reg_blockMatching_gpu)

# Build NiftySeg libraries: 
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/seg-lib)
SUBDIRS(seg-lib)

# Build Matlab extensions
IF(NiftyRec_BUILD_MATLAB)
        SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake/matlab")
        FIND_PACKAGE(MATLAB REQUIRED)
        SUBDIRS(Matlab)
ENDIF(NiftyRec_BUILD_MATLAB)

# Build volume renderer
IF(NiftyRec_BUILD_RENDERER)
        SUBDIRS(volume_render)
ENDIF(NiftyRec_BUILD_RENDERER)

# Build documentation
IF(NiftyRec_COMPILE_DOC)
        SUBDIRS(documentation)
ENDIF(NiftyRec_COMPILE_DOC)

# Build applications
IF(NiftyRec_BUILD_APPS)
        SUBDIRS(apps)
ENDIF(NiftyRec_BUILD_APPS)

# Build data 
IF(NiftyRec_INCLUDE_DATA)
    SUBDIRS(data)
ENDIF(NiftyRec_INCLUDE_DATA)

# Build miscellaneous 
SUBDIRS(misc)

#-----------------------------------------------------------------------------
# Set installation directories: 
SET(LIB_INSTALL_DIR "lib")
SET(BIN_INSTALL_DIR "niftyrec/bin")
SET(INC_INSTALL_DIR "include")
SET(DOC_INSTALL_DIR "niftyrec/documentation")
SET(MEX_INSTALL_DIR_NIFTYREC "niftyrec/matlab")
SET(MEX_INSTALL_DIR_NIFTYREG "niftyreg/matlab")
SET(MEX_INSTALL_DIR_NIFTYSEG "niftyseg/matlab")

#-----------------------------------------------------------------------------
# Build NiftyRec packages: 
# 1) .deb installer for Linux Debian, Ubuntu, Knoppix, ..; 
# 2) NSIS self installer for Windows;
# 3) Drag and drop installer for MacOsX;
# 4) Compressed source tree. 

IF(NiftyRec_DO_PACKAGE)
        IF(UNIX)
            IF(APPLE)
                SET(CPACK_GENERATOR "TGZ;PackageMaker")
            ELSE(APPLE)
                SET(CPACK_GENERATOR "TGZ;DEB")
            ENDIF(APPLE)
        ENDIF(UNIX)

        SET(CPACK_PACKAGE_NAME "NiftyRec")
        SET(CPACK_PACKAGE_VENDOR "CMIC-UCL")
        SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "NiftyRec - Tomographic Reconstruction")
        SET(CPACK_PACKAGE_VERSION_MAJOR ${Nifty_Rec_VERSION_MAJOR})
        SET(CPACK_PACKAGE_VERSION_MINOR ${Nifty_Rec_VERSION_MINOR})
        SET(CPACK_PACKAGE_VERSION_PATCH ${Nifty_Rec_VERSION_PATCH})
        SET(CPACK_PACKAGE_INSTALL_DIRECTORY "NiftyRec")
        SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.txt")
        SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.txt")
        SET(CPACK_PACKAGE_ICON "${CMAKE_CURRENT_SOURCE_DIR}/images/niftyrec_logo.ico")
        SET(CPACK_NSIS_MUI_ICON "${CMAKE_CURRENT_SOURCE_DIR}/images/niftyrec_logo.ico")
        SET(CPACK_NSIS_MUI_UNIICON "${CMAKE_CURRENT_SOURCE_DIR}/images/niftyrec_logo.ico")

        IF (CMAKE_SYSTEM_PROCESSOR MATCHES "unknown")
          SET (CMAKE_SYSTEM_PROCESSOR "x86")
        ENDIF (CMAKE_SYSTEM_PROCESSOR MATCHES "unknown")
        IF(NOT DEFINED CPACK_SYSTEM_NAME)
          SET(CPACK_SYSTEM_NAME ${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
        ENDIF(NOT DEFINED CPACK_SYSTEM_NAME)
        IF(${CPACK_SYSTEM_NAME} MATCHES Windows)
          IF(CMAKE_CL_64)
            SET(CPACK_SYSTEM_NAME Win64-${CMAKE_SYSTEM_PROCESSOR})
          ELSE(CMAKE_CL_64)
            SET(CPACK_SYSTEM_NAME Win32-${CMAKE_SYSTEM_PROCESSOR})
          ENDIF(CMAKE_CL_64)
        ENDIF(${CPACK_SYSTEM_NAME} MATCHES Windows)

        IF(${CPACK_SYSTEM_NAME} MATCHES Darwin AND CMAKE_OSX_ARCHITECTURES)
          list(LENGTH CMAKE_OSX_ARCHITECTURES _length)
          IF(_length GREATER 1)
            SET(CPACK_SYSTEM_NAME Darwin-Universal)
          ELSE(_length GREATER 1)
            SET(CPACK_SYSTEM_NAME Darwin-${CMAKE_OSX_ARCHITECTURES})
          ENDIF(_length GREATER 1)
        ENDIF(${CPACK_SYSTEM_NAME} MATCHES Darwin AND CMAKE_OSX_ARCHITECTURES)

        # .deb package
        SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "Stefano Pedemonte stefano.pedemonte@gmail.com")

        # Nsis package
        SET(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\NiftyRec.exe")
        SET(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY} NiftyRec")
        SET(CPACK_NSIS_HELP_LINK "http:\\\\\\\\niftyrec.scienceontheweb.net")
        SET(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\niftyrec.scienceontheweb.net")
        SET(CPACK_NSIS_CONTACT "stefano.pedemonte@gmail.com")
        SET(CPACK_NSIS_MODIFY_PATH ON)

        # Apple drag-and-drop package
        IF(APPLE)
                SET(CPACK_PACKAGE_INSTALL_DIRECTORY "/Library/niftyrec/")
                SET(CPACK_PACKAGING_INSTALL_PREFIX. "/Library/niftyrec/")
                SET(CPACK_BINARY_TBZ2 OFF)
                SET(CPACK_BINARY_DRAGNDROP ON)
                SET(CPACK_BINARY_PACKAGEMAKER ON)
                SET(CPACK_BINARY_STGZ OFF)
        ENDIF(APPLE)

        INCLUDE(CPack)
ENDIF(NiftyRec_DO_PACKAGE)

#-----------------------------------------------------------------------------
# Hide unnecessary build options 
MACRO ( MARK_AS_INTERNAL _var )
  SET ( ${_var} ${${_var}} CACHE INTERNAL "hide this!" FORCE )
ENDMACRO( MARK_AS_INTERNAL _var ) 

MARK_AS_INTERNAL(Teem_PNG)
MARK_AS_INTERNAL(Teem_PTHREAD)
MARK_AS_INTERNAL(Teem_USE_LIB_INSTALL_SUBDIR)
MARK_AS_INTERNAL(Teem_VTK_MANGLE)
MARK_AS_INTERNAL(Teem_VTK_TOOLKITS_IPATH)
MARK_AS_INTERNAL(Teem_VTK_ZLIB_MANGLE_IPATH)
MARK_AS_INTERNAL(Teem_ZLIB)
MARK_AS_INTERNAL(Teem_BZIP2)
MARK_AS_INTERNAL(Teem_FFTW3)
MARK_AS_INTERNAL(Teem_LEVMAR)
MARK_AS_INTERNAL(BUILD_SHARED_LIBS)
MARK_AS_INTERNAL(BUILD_HEX)
MARK_AS_INTERNAL(BUILD_EXPERIMENTAL_APPS)
MARK_AS_INTERNAL(BUILD_EXPERIMENTAL_LIBS)
MARK_AS_INTERNAL(CUDA_BUILD_CUBIN) 
MARK_AS_INTERNAL(CUDA_VERBOSE_BUILD) 
MARK_AS_INTERNAL(CUDA_SDK_ROOT_DIR) 
#----------------------------------------------------------------------------- 


