lptModule(dataaquisition core calib datagen imageproc correspond tracking visualization)

set (EXTRA_LIBS ${EXTRA_LIBS} ${OpenCV_LIBS} ${VTK_LIBRARIES} )
if (NP_CAMERASDK_FOUND)
	add_definitions(-DUSE_NP_CAMERASDK)
	set (EXTRA_LIBS ${EXTRA_LIBS} ${NP_CAMERASDK_LIBRARIES} )
endif()

if (CUDA_FOUND)
	include_directories("${CMAKE_SOURCE_DIR}/modules/gpu/include")
	target_link_libraries(${PROJECT_NAME} module-correspondcuda)
endif()

target_link_libraries(${PROJECT_NAME} ${EXTRA_LIBS} )