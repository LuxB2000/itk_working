# This project is designed to be built outside the Insight source tree.
PROJECT(itkSuarezBlockMatchingRegistration2)

INCLUDE(${CMAKE_SOURCE_DIR}/IJMacros.txt)

ENABLE_TESTING()
INCLUDE(Dart)

SET(Required_Packages
ITK2.8.1
)

FOREACH(Package ${Required_Packages})
  LOADPACKAGE(${Package})
ENDFOREACH(Package)

SET(Libraries
ITKIO
ITKBasicFilters
ITKNumerics
ITKCommon
ITKAlgorithms
ITKStatistics
)

SET(CurrentExe "ImageCompare")
ADD_EXECUTABLE(${CurrentExe} ImageCompare.cxx) 
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})

SET(CurrentExe "itkSuarezBlockMatchingRegistration2")
ADD_EXECUTABLE(${CurrentExe} main.cxx)
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})

ADD_TEST(RegisterImage itkSuarezBlockMatchingRegistration2 ${CMAKE_SOURCE_DIR}/data_input0_fevbsep256.jpg  ${CMAKE_SOURCE_DIR}/data_input1_odraude256.jpg 2.0 3.0 0 ${CMAKE_SOURCE_DIR}/data_outputMine_odraude256def.mha)
ADD_TEST(CompareImage ImageCompare ${CMAKE_SOURCE_DIR}/data_outputAuthor_odraude256def.mha ${CMAKE_SOURCE_DIR}/data_outputMine_odraude256def.mha)

