project("module-correspondcuda")
message(STATUS "Adding module correspondcuda")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include"
                        "${CMAKE_CURRENT_SOURCE_DIR}/src"
                        "${CMAKE_CURRENT_BINARY_DIR}")
        
include_directories("${CMAKE_SOURCE_DIR}/modules/core/include")
include_directories("${CMAKE_SOURCE_DIR}/modules/correspond/include")

message(STATUS "   depends on module: core")
message(STATUS "   depends on module: correspond")

file(GLOB lib_int_hdrs "src/*.h*")
file(GLOB lib_cu_srcs "src/*.cu")
file(GLOB lib_cpp_srcs "src/*.cpp")
source_group("Src" FILES ${lib_cpp_srcs} ${lib_cu_srcs} ${lib_int_hdrs})


file(GLOB lib_hdrs "include/*.h*")
source_group("Include" FILES ${lib_hdrs})
 
cuda_add_library(${PROJECT_NAME} ${lib_cu_srcs} ${lib_cpp_srcs} ${lib_hdrs} ${lib_int_hdrs})
    
target_link_libraries(${PROJECT_NAME} module-core module-correspond ${CUDA_LIBRARIES})
add_dependencies(${PROJECT_NAME} module-core module-correspond)

if(Boost_FOUND)
	target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})
endif()   
