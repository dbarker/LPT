#Lagrangian particle tracking module definition

macro(lptModule name)
    
    project("module-${name}")
    message(STATUS "Adding module ${name}")
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include"
                        "${CMAKE_CURRENT_SOURCE_DIR}/src"
                        "${CMAKE_CURRENT_BINARY_DIR}")
        
    
    foreach(lpt_module ${ARGN})
        include_directories("${CMAKE_SOURCE_DIR}/modules/${lpt_module}/include")
    endforeach()

    file(GLOB lib_srcs "src/*.cpp")
    file(GLOB lib_int_hdrs "src/*.h*")
    source_group("Src" FILES ${lib_srcs} ${lib_int_hdrs})

    file(GLOB lib_hdrs "include/*.h*")
    source_group("Include" FILES ${lib_hdrs})
 
    add_library(${PROJECT_NAME} ${lib_srcs} ${lib_hdrs} ${lib_int_hdrs})
    
    foreach(lpt_module ${ARGN})
        set(module_dep "module-${lpt_module}")
        target_link_libraries(${PROJECT_NAME} ${module_dep})
        add_dependencies(${PROJECT_NAME} ${module_dep})
        message(STATUS "   depends on module: ${lpt_module}")
    endforeach()
       
    if(Boost_FOUND)
	target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})
    endif()   
            
endmacro()

macro(lptApplication name)
  project("app-${name}")
endmacro()

macro(lptTest name)
  project("test-${name}")
endmacro()

