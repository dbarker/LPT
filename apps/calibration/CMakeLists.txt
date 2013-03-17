#CMakeLists file for samples/calib
message(STATUS "generating makefiles for calibration: ${CMAKE_SOURCE_DIR}")
lptApplication(Calibration)
include_directories("${CMAKE_SOURCE_DIR}/modules/core/include")
include_directories("${CMAKE_SOURCE_DIR}/modules/calib/include")

set (EXTRA_LIBS ${EXTRA_LIBS} module-core module-calib)

add_executable(${PROJECT_NAME} calibrator.cpp)
target_link_libraries(${PROJECT_NAME} ${EXTRA_LIBS} )