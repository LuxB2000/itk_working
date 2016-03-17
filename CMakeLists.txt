cmake_minimum_required(VERSION 2.6)
 
PROJECT(IETools)

# Specify the path if you use a non default ITK location
# SET(ITK_DIR "/mine/InsightToolkit-4.8.1-build")

FIND_PACKAGE ( ITK REQUIRED )
IF (  ITK_FOUND  )
	INCLUDE(  ${ITK_USE_FILE}  )
	MESSAGE( ">> ITK files: ${ITK_USE_FILE} ")
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR "ITK not found. Please set ITK_DIR")
ENDIF(   ITK_FOUND   )

# Specify the path if you use a non default ITK location
# SET(ITK_DIR "/mine/VTK")

FIND_PACKAGE ( VTK )
IF (  VTK_FOUND  )
	INCLUDE(  ${VTK_USE_FILE}  )
	MESSAGE( ">> VTK files: ${VTK_USE_FILE} ")
ELSE( VTK_FOUND )
  MESSAGE(FATAL_ERROR "VTK not found. Please set VTK_DIR")
ENDIF(   VTK_FOUND   )

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# TARGET_LINK_LIBRARIES( IETools ITKCommon ITKIO )

add_subdirectory( tools )
