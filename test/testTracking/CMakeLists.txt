message(STATUS "generating makefiles for TestTracking")
lptTest(Tracking)
include_directories("${CMAKE_SOURCE_DIR}/modules/core/include")
include_directories("${CMAKE_SOURCE_DIR}/modules/tracking/include")

set (EXTRA_LIBS ${EXTRA_LIBS} module-core module-tracking)

add_executable(${PROJECT_NAME} trackingtest.cpp)
target_link_libraries(${PROJECT_NAME} ${EXTRA_LIBS} )