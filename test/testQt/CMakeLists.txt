message(STATUS "generating makefiles for qt test")
lptTest(Qt)

include_directories( ${Boost_INCLUDE_DIRS} )
include_directories("${CMAKE_SOURCE_DIR}/modules/core/include")

set (EXTRA_LIBS ${EXTRA_LIBS} module-core ${Boost_LIBRARIES})

set(qttest_SOURCES qttest.cpp )
set(qttest_HEADERS )

QT4_WRAP_CPP(qttest_HEADERS_MOC ${qttest_HEADERS})
include(${QT_USE_FILE})
add_definitions(${QT_DEFINITIONS})

add_executable(${PROJECT_NAME}  ${qttest_SOURCES} ${qttest_HEADERS_MOC})
target_link_libraries(${PROJECT_NAME}  ${EXTRA_LIBS} )
target_link_libraries(${PROJECT_NAME}  ${QT_LIBRARIES})