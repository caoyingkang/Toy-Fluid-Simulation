//
// Created by cyk on 18-12-27.
//

#ifndef PROJ_SHADER_H
#define PROJ_SHADER_H

#include <glad/glad.h> // include glad to get all the required OpenGL headers
#include <string>

class Shader {
public:
    // constructor reads and builds the shader
    Shader(const GLchar *vertexPath, const GLchar *fragmentPath);

    // use/activate the shader
    void use() {
        glUseProgram(ID);
    }

    // utility uniform functions
    void setBool(const std::string &name, bool value) const {
        glUniform1i(glGetUniformLocation(ID, name.c_str()), (int) value);
    }

    void setInt(const std::string &name, int value) const {
        glUniform1i(glGetUniformLocation(ID, name.c_str()), value);
    }

    void setFloat(const std::string &name, float value) const {
        glUniform1f(glGetUniformLocation(ID, name.c_str()), value);
    }

private:
    // the program ID
    unsigned int ID;
};


#endif //PROJ_SHADER_H
