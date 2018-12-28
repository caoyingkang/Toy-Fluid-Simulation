#version 330 core

in vec4 currColor; // r,b
out vec4 FragColor;

void main()
{
    FragColor = currColor;
}