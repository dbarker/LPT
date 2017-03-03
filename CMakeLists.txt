# Main CMakeLists file

cmake_minimum_required (VERSION 2.8)
project (LPT)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/scripts/cmake/")
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
  
include( "./scripts/cmake/LptCmakeMacros.cmake" REQUIRED)  
find_package( OpenCV REQUIRED )
find_package( VTK 5.10 REQUIRED)
find_package( YamlCPP REQUIRED )
find_package( Qt4 )
set( QT_USE_QTOPENGL TRUE )

find_package( CUDA 4.2 )

message(STATUS "OpenCV include: " ${OpenCV_INCLUDE_DIRS})

link_directories( ${VTK_LIBRARY_DIRS} ) 
include_directories( ${VTK_INCLUDE_DIRS} )
message(STATUS "VTK include: " ${VTK_INCLUDE_DIRS})
message(STATUS "VTK lib: " ${VTK_LIBRARY_DIRS})

if (CUDA_FOUND)
	link_directories( ${CUDA_LIBRARY_DIRS} )
	include_directories( ${CUDA_INCLUDE_DIRS} )
	add_definitions( -DUSE_CUDA )
	message(STATUS "CUDA include: " ${CUDA_INCLUDE_DIRS})
	message(STATUS "CUDA toolkit: " ${CUDA_TOOLKIT_ROOT_DIR})
endif()

set(Boost_USE_STATIC_LIBS        ON)
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME     OFF)
set(Boost_NO_BOOST_CMAKE         ON)
set(Boost_ADDITIONAL_VERSIONS "1.48" "1.48.0")
set(Boost_INCLUDE_DIR $ENV{Boost_INCLUDE_DIR})

find_package( Boost COMPONENTS thread chrono random date_time system REQUIRED )
message(STATUS "Boost information: BOOST_ROOT = $ENV{BOOST_ROOT}") 
message(STATUS "Boost_INCLUDE_DIRS: ${Boost_INCLUDE_DIRS}") 
message(STATUS "Boost_LIBRARY_DIRS: ${Boost_LIBRARY_DIRS}") 

if (WIN32)
    find_package( NP_CAMERASDK )
	#add_definitions(${Boost_LIB_DIAGNOSTIC_DEFINITIONS})
    include_directories($ENV{MSMPI_INC})
endif()

include_directories(${Boost_INCLUDE_DIRS} ${YAMLCPP_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIRS})

message(STATUS "Starting LPT library build:")
add_subdirectory("${CMAKE_SOURCE_DIR}/modules")
add_subdirectory("${CMAKE_SOURCE_DIR}/apps")
add_subdirectory("${CMAKE_SOURCE_DIR}/test")