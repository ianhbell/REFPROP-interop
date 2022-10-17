project(interop)

set(CMAKE_CXX_STANDARD 17)

add_library(RPinterop INTERFACE)
target_include_directories(RPinterop INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/externals/args_hxx")
target_include_directories(RPinterop INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include")

add_executable(main "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp")
target_link_libraries(main PRIVATE RPinterop)