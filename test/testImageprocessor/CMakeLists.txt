lptTest(ImageProcessor)
include_directories("${CMAKE_SOURCE_DIR}/modules/core/include")
include_directories("${CMAKE_SOURCE_DIR}/modules/imageproc/include")
include_directories("${CMAKE_SOURCE_DIR}/modules/datagen/include")
include_directories("${CMAKE_SOURCE_DIR}/modules/correspond/include")

set (EXTRA_LIBS 
    module-core 
	module-imageproc 
	module-datagen
	)

add_executable(${PROJECT_NAME}  imageprocessortest.cpp)
target_link_libraries( ${PROJECT_NAME}  ${EXTRA_LIBS} )