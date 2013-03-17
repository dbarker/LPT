#CMakeLists file for concurrent queue test

lptTest(ConcurrentQueue)
include_directories("${CMAKE_SOURCE_DIR}/modules/core/include")

set (EXTRA_LIBS ${EXTRA_LIBS} module-core ${OpenCV_LIBS})

add_executable(${PROJECT_NAME}  concurrentQueueTest.cpp)
target_link_libraries( ${PROJECT_NAME}  ${EXTRA_LIBS} )