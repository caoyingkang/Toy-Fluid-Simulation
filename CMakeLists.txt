cmake_minimum_required(VERSION 3.0)
project(proj)

set(CMAKE_CXX_STANDARD 11)
set(GLADCPATH /home/cyk/openglDeps/glad/src/glad.c)

add_executable(proj
        main.cpp Simulator.cpp
        Simulator.h dCmpUtil.h
        ${GLADCPATH})

target_link_libraries(proj
        GLEW glfw3 dl GL X11
        pthread Xrandr Xxf86vm
        Xinerama Xcursor Xi)

add_executable(test main2.cpp)