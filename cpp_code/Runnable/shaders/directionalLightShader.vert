#version 330 core
// NOTE: Do NOT use any version older than 330! Bad things will happen!

// This is an example vertex shader. GLSL is very similar to C.
// You can define extra functions if needed, and the main() function is
// called when the vertex shader gets run.
// The vertex shader gets called once per vertex.

layout (location = 0) in vec3 position;
layout (location = 1) in vec3 normal;

// Uniform variables can be updated by fetching their location and passing values to that location
uniform mat4 projection;
uniform mat4 view;
uniform mat4 model;
uniform vec3 eye;

// Matrices used to transform normal variables
uniform mat4 projectionTI;
uniform mat4 viewTI;
uniform mat4 modelTI;

// Information from the material
uniform vec3 color;
uniform vec3 ambient;
uniform vec3 diffuse;
uniform vec3 specular;
uniform float shininess;

// Information from the light
uniform vec3 lightDir;
uniform mat3 lightColor;

// Outputs of the vertex shader are the inputs of the same name of the fragment shader.
// The default output, gl_Position, should be assigned something. You can define as many
// extra outputs as you need.
out vec3  colorInfo;
out float alphaValue;

void main()
{
    // OpenGL maintains the D matrix so you only need to multiply by P, V (aka C inverse), and M
    vec4 vertexPosition4 = model * vec4(position, 1.0);
    gl_Position = projection * view * vertexPosition4;
    
    vec3 vertexPosition = vec3(vertexPosition4[0], vertexPosition4[1], vertexPosition4[2]);
    
    // Calculate color with light
    vec3 lightDirNor = normalize(lightDir);
    vec4 normal4     = modelTI * vec4(normal, 1.0);
    vec3 normalNor   = normalize(vec3(normal4[0], normal4[1], normal4[2]));
    float diffuseCo = max(dot(normalNor, lightDirNor), 0.0);

    // Specular: -gl_Position is the viewing directory (?)
    vec3 reflectDir = normalize(reflect(lightDirNor, normalNor));
    vec3 viewDir = normalize(vertexPosition-eye);
    float specularCo = pow(max(dot(viewDir, reflectDir), 0.0), shininess);

    colorInfo = lightColor * ambient + diffuseCo * lightColor * diffuse + specularCo * lightColor * specular;
    alphaValue  = 1.0f;
}