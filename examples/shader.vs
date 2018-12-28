#version 330 core

layout (location=0) in vec2 aPos; // x,y
layout (location=1) in vec2 aColor; // r,b

out vec4 currColor;

void main()
{
    gl_Position = vec4(aPos, 0.0, 1.0);
    currColor = vec4(aColor.x, 0.0, aColor.y, 1.0);
}
