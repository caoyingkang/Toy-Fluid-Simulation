cmake_minimum_required(VERSION 3.0)
project(ToyFluidSimProj)

set(CMAKE_CXX_STANDARD 11)

if (CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "-O2")        ## Optimize
    set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary
endif ()

set(GLADCPATH /home/cyk/openglDeps/glad/src/glad.c)

configure_file(shader.vs . COPYONLY)
configure_file(shader.fs . COPYONLY)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(ToyFluidSimProj
        main.cpp
        cmpUtil.h
        Simulator.h Simulator.cpp
        Shader.h Shader.cpp
        cgUtil.h cgUtil.cpp
        ${GLADCPATH})

target_link_libraries(ToyFluidSimProj
        GLEW glfw3 dl GL X11
        pthread Xrandr Xxf86vm
        Xinerama Xcursor Xi
        Eigen3::Eigen)


