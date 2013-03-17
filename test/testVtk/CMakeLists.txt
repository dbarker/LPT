message(STATUS "generating makefiles for vtk test")
lptTest(visualization)

include_directories( ${Boost_INCLUDE_DIRS} )
include_directories( ${VTK_INCLUDE_DIRS} )
include_directories("${CMAKE_SOURCE_DIR}/modules/core/include")
include_directories("${CMAKE_SOURCE_DIR}/modules/visualization/include")

set (EXTRA_LIBS ${EXTRA_LIBS} module-core module-visualization ${VTK_LIBRARIES} ${Boost_LIBRARIES})

add_executable(${PROJECT_NAME} vtktest.cpp)
target_link_libraries(${PROJECT_NAME} ${EXTRA_LIBS} )