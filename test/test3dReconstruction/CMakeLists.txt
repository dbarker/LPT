lptTest(3dReconstruction)

include_directories( ${Boost_INCLUDE_DIRS} )
include_directories("${CMAKE_SOURCE_DIR}/modules/core/include")
include_directories("${CMAKE_SOURCE_DIR}/modules/correspond/include")
include_directories("${CMAKE_SOURCE_DIR}/modules/datagen/include")

set (EXTRA_LIBS ${EXTRA_LIBS} 
	module-core 
	module-correspond 
	module-datagen 
	)

if(CUDA_FOUND)
	include_directories("${CMAKE_SOURCE_DIR}/modules/gpu/include")
	set(EXTRA_LIBS ${EXTRA_LIBS} module-correspondcuda )
endif()

add_executable(${PROJECT_NAME} reconstruct3dtest.cpp)
target_link_libraries(${PROJECT_NAME} ${EXTRA_LIBS} )