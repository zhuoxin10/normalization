
PROJECT(Normalization)

cmake_minimum_required(VERSION 2.6)
CMAKE_POLICY(SET CMP0015 NEW)

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
          "ITK not found. Please set ITK_DIR.")
ENDIF(ITK_FOUND)

include_directories(${CMAKE_CURRENT_BINARY_DIR})
link_directories(${CMAKE_CURRENT_BINARY_DIR})

ADD_DEFINITIONS(-DNDEBUG -Wall -g -c -std=gnu++11)

ADD_EXECUTABLE(normalization normalization.cc)

TARGET_LINK_LIBRARIES(normalization ${ITK_LIBRARIES} )

