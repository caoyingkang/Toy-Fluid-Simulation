//
// Created by cyk on 18-12-23.
//


#include "Simulator.h"
#include "Shader.h"
#include <iostream>
#include <iomanip>
#include <glad/glad.h>
#include <GLFW/glfw3.h>

using namespace std;

// settings for simulation
const float dx = 0.1;
const float dt = 0.1;
const int Nx = 100;
const int Ny = 100;
const float visc = 0;
const float diff = 0;

// settings for rendering
const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 600;
const char vsFile[] = "shader.vs";
const char fsFile[] = "shader.fs";
const int n_triangles = 2 * (Nx - 1) * (Ny - 1);
const int sz_vertex = 4; // x,y,(z=0),r,(g=0),b

// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------
void processInput(GLFWwindow *window);

// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
void framebuffer_size_callback(GLFWwindow *window, int width, int height);


void setPosVertices(float *vertices);


int main() {
    // initializing simulator
    // ----------------------
    Simulator smlt(dx, dt, {Nx, Ny}, visc, diff);

    // set force
    // default: no force
    // smlt.setForce(0,10);
    // smlt.setForce(some VecMatXd);

    // set inlet
    // smlt.setDefaultInLet();
    // smlt.setInlet(radius_blue, {center_blue}, {v_blue},
    //               radius_red, {center_red}, {v_red})
//    smlt.setInlet(6, {25, 25}, {10, 10},
//                  6, {75, 25}, {-10, 10});
    smlt.setInlet(3, {25, 25}, {1, 1},
                  3, {25, 75}, {1.1, -1.1});

    // glfw: initialize and configure
    // ------------------------------
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // uncomment this statement to fix compilation on OS X
#endif

    // glfw window creation
    // --------------------
    GLFWwindow *window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "Toy-Liquid-Simulation", nullptr, nullptr);
    if (window == nullptr) {
        cerr << "Failed to create GLFW window" << endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    // glad: load all OpenGL function pointers
    // ---------------------------------------
    if (!gladLoadGLLoader((GLADloadproc) glfwGetProcAddress)) {
        cerr << "Failed to initialize GLAD" << endl;
        return -1;
    }

    // build and compile our shader program
    // ------------------------------------
    Shader shd(vsFile, fsFile);

    // set up vertex data (and buffer(s)) and configure position attributes
    // ------------------------------------------------------------------
    float vertices[3 * n_triangles * sz_vertex];
    setPosVertices(vertices);

    unsigned int VBO, VAO;
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);

    // bind the Vertex Array Object first,
    // then bind and set vertex buffer(s),
    // and then configure vertex attributes(s).
    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_DYNAMIC_DRAW);

    // position attribute (location=0)
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void *) 0);
    glEnableVertexAttribArray(0);
    // color attribute (location=1)
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void *) (2 * sizeof(float)));
    glEnableVertexAttribArray(1);

    // render loop
    // -----------
    int MaxIter = 300;
    int iter = 0;
    while (!glfwWindowShouldClose(window)) {
        // input
        // -----
        processInput(window);

        // simulation step
        if (iter < MaxIter) {
            smlt.Forward();
            smlt.getRenderData(vertices);
            //
            smlt.Forward();
            //
            if (iter + 1 == MaxIter) {
                smlt.printBlue();
                smlt.printRed();
                smlt.printVx();
                smlt.printVy();
                int idx = 0;
                cout << "vertices color:" << endl;
                while (idx < 3 * n_triangles * sz_vertex) {
                    cout << setw(16) << setprecision(4) << fixed
                         << vertices[idx] << "\t\t\t"
                         << vertices[idx + 1] << "\t\t\t"
                         << vertices[idx + 2] << "\t\t\t"
                         << vertices[idx + 3] << endl;
                    idx += 4;
                }
            }
            //
            ++iter;
        }


        // render
        // ------
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        // render the triangle
        shd.use();
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices);
        glBindVertexArray(VAO);
        glDrawArrays(GL_TRIANGLES, 0, 3 * n_triangles);

        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        // -------------------------------------------------------------------------------
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // optional: de-allocate all resources once they've outlived their purpose:
    // ------------------------------------------------------------------------
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);

    // glfw: terminate, clearing all previously allocated GLFW resources.
    // ------------------------------------------------------------------
    glfwTerminate();
    return 0;
}


void processInput(GLFWwindow *window) {
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}

void framebuffer_size_callback(GLFWwindow *window, int width, int height) {
    // make sure the viewport matches the new window dimensions; note that width and
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
}

void setPosVertices(float *vertices) {
    int i, j;
    int idx = 0; // index of vertices
    float deltax = 2.0f / Nx, deltay = 2.0f / Ny;
    for (j = 0; j < Ny - 1; ++j) {
        for (i = 0; i < Nx - 1; ++i) {
            float posx = 2 * float(i) / Nx - 1,
                    posy = 2 * float(j) / Ny - 1;
            vertices[idx] = posx;
            vertices[idx + 1] = posy;
            idx += 4;
            vertices[idx] = posx + deltax;
            vertices[idx + 1] = posy;
            idx += 4;
            vertices[idx] = posx;
            vertices[idx + 1] = posy + deltay;
            idx += 4;
            vertices[idx] = posx + deltax;
            vertices[idx + 1] = posy;
            idx += 4;
            vertices[idx] = posx;
            vertices[idx + 1] = posy + deltay;
            idx += 4;
            vertices[idx] = posx + deltax;
            vertices[idx + 1] = posy + deltay;
            idx += 4;
        }
    }
}