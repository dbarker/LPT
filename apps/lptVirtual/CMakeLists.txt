lptApplication(LPT_Virtual)

message(STATUS "generating makefiles ${PROJECT_NAME}")

include_directories( ${Boost_INCLUDE_DIRS} )
include_directories("${CMAKE_SOURCE_DIR}/modules/core/include")
include_directories("${CMAKE_SOURCE_DIR}/modules/dataaquisition/include")
include_directories("${CMAKE_SOURCE_DIR}/modules/imageproc/include")
include_directories("${CMAKE_SOURCE_DIR}/modules/datagen/include")
include_directories("${CMAKE_SOURCE_DIR}/modules/calib/include")
include_directories("${CMAKE_SOURCE_DIR}/modules/correspond/include")
include_directories("${CMAKE_SOURCE_DIR}/modules/tracking/include")
include_directories("${CMAKE_SOURCE_DIR}/modules/visualization/include")

set (EXTRA_LIBS 
    module-core
	module-calib
    module-imageproc 
    module-correspond 
    module-tracking 
    module-datagen 
    module-dataaquisition 
    module-visualization 
    )

if (CUDA_FOUND)
    include_directories("${CMAKE_SOURCE_DIR}/modules/gpu/include")
    set(EXTRA_LIBS ${EXTRA_LIBS} module-correspondcuda )
endif()

add_executable(${PROJECT_NAME} virtualcams.cpp)
target_link_libraries(${PROJECT_NAME} ${EXTRA_LIBS} )