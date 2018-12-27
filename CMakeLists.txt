cmake_minimum_required(VERSION 3.0)
project(proj)

set(CMAKE_CXX_STANDARD 11)
set(GLADCPATH /home/cyk/openglDeps/glad/src/glad.c)

FIND_PACKAGE(OpenMP REQUIRED)
if (OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif ()


add_executable(proj
        main.cpp
        dCmpUtil.h
        Simulator.h Simulator.cpp
        Shader.h Shader.cpp
        cgUtil.h cgUtil.cpp
        ${GLADCPATH})

target_link_libraries(proj
        GLEW glfw3 dl GL X11
        pthread Xrandr Xxf86vm
        Xinerama Xcursor Xi)

add_executable(test main2.cpp)