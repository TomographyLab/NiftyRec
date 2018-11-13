

find_package(CUDA)
find_package(GLUT)
find_package(OpenGL)

OPTION(FORCE_32_BIT      "set 32 bit copiler options" OFF)
SET(VOLUMERENDER_INSTALL_DIR "volume_render")

IF (WIN32)
    FIND_PATH( GLEW_INCLUDE_PATH GL/glew.h
        $ENV{PROGRAMFILES}/GLEW/include
        ${PROJECT_SOURCE_DIR}/src/nvgl/glew/include
        DOC "The directory where GL/glew.h resides")
    FIND_LIBRARY( GLEW_LIBRARY
        NAMES glew GLEW glew32 glew32s
        PATHS
        $ENV{PROGRAMFILES}/GLEW/lib
        ${PROJECT_SOURCE_DIR}/src/nvgl/glew/bin
        ${PROJECT_SOURCE_DIR}/src/nvgl/glew/lib
        DOC "The GLEW library")
ELSE (WIN32)
    FIND_PATH( GLEW_INCLUDE_PATH GL/glew.h
        /usr/include
        /usr/local/include
        /sw/include
        /opt/local/include
        DOC "The directory where GL/glew.h resides")
    FIND_LIBRARY( GLEW_LIBRARY
        NAMES GLEW glew
        PATHS
        /usr/lib64
        /usr/local/lib64
        DOC "The GLEW library")
ENDIF (WIN32)

IF (GLEW_INCLUDE_PATH)
    SET( GLEW_FOUND 1 CACHE STRING "Set to 1 if GLEW is found, 0 otherwise")
ELSE (GLEW_INCLUDE_PATH)
    SET( GLEW_FOUND 0 CACHE STRING "Set to 1 if GLEW is found, 0 otherwise")
ENDIF (GLEW_INCLUDE_PATH)
MESSAGE(STATUS,"Glew library: ${GLEW_LIBRARY}")
MESSAGE(STATUS,"Glew include: ${GLEW_INCLUDE_PATH}")
MARK_AS_ADVANCED( GLEW_FOUND )

#FIND_PATH( SHRUTIL_INCLUDE_PATH shrUtils.h
#        ${CUDA_SDK_ROOT_DIR}/../shared/inc
#        DOC "The directory where shrutil.h resides")
#FIND_LIBRARY( SHRUTIL_LIBRARY
#        NAMES SHRUTIL shrutil_i386_64
#        PATHS
#        ${CUDA_SDK_ROOT_DIR}/../shared/lib
#        DOC "The SHRUTIL library")

#SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE) 
#SET(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX})SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE) 
#SET(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX})

include(FindCUDA)
set(CORELIBS 
    ${GLUT_LIBRARY} 
    ${OPENGL_LIBRARY}
    ${GLEW_LIBRARY} 
    ${CUDA_SDK_LIBRARIES} 
    m)

IF(FORCE_32_BIT)
    set(CUDA_64_BIT_DEVICE_CODE OFF) 
    set(CMAKE_CXX_FLAGS "-g -Wall -m32")
    if (Apple)
        set (CMAKE_OSX_ARCHITECTURES i386)
    endif (Apple)
ENDIF(FORCE_32_BIT)

link_libraries(${CORELIBS})

include_directories(
    ${CUDA_SDK_ROOT_DIR}/common/inc
    ${CUDA_SDK_ROOT_DIR}/../shared/inc
    ${CMAKE_CURRENT_SOURCE_DIR}
)

CUDA_ADD_LIBRARY(volumeRender SHARED volumeRender.h volumeRender.cpp volumeRender_kernel.h volumeRender_kernel.cu)
INSTALL(TARGETS volumeRender DESTINATION ${VOLUMERENDER_INSTALL_DIR})
#CUDA_ADD_EXECUTABLE(volumeRenderTest volumeRender.h volumeRenderTest.cpp )
#INSTALL(TARGETS volumeRenderBin DESTINATION ${VOLUMERENDER_INSTALL_DIR})
#TARGET_LINK_LIBRARIES(volumeRenderBin volumeRender)

#install Python files
FILE(GLOB PY_FILES "*.py")
INSTALL(FILES ${PY_FILES} DESTINATION ${VOLUMERENDER_INSTALL_DIR})
#install shell scripts
FILE(GLOB SH_FILES "*.sh")
INSTALL(FILES ${SH_FILES} DESTINATION ${VOLUMERENDER_INSTALL_DIR})
#install example .mat files
FILE(GLOB MAT_FILES "*.mat")
INSTALL(FILES ${MAT_FILES} DESTINATION ${VOLUMERENDER_INSTALL_DIR})



