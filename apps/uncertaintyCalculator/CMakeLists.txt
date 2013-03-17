#CMakeLists file for uncertainty calculator

lptApplication(UncertaintyCalculator)
include_directories("${CMAKE_SOURCE_DIR}/modules/core/include")
include_directories("${CMAKE_SOURCE_DIR}/modules/calib/include")
include_directories("${CMAKE_SOURCE_DIR}/modules/datagen/include")

set (EXTRA_LIBS ${EXTRA_LIBS}
 module-core 
 module-datagen 
 module-calib
 )

add_executable(${PROJECT_NAME} calculateuncertainties.cpp)
target_link_libraries( ${PROJECT_NAME} ${EXTRA_LIBS} )