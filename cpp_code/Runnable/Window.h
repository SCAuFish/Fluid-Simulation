#ifndef _WINDOW_H_
#define _WINDOW_H_

#ifdef __APPLE__
#define GLFW_INCLUDE_GLCOREARB
#include <OpenGL/gl3.h>
#else
#include <GL/glew.h>
#endif
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <iostream>
#include <vector>
#include <string>
#include <memory>

#include "Object.h"
#include "Cube.h"
#include "PointCloud.h"
#include "shader.h"
#include "RasterizerQuad.h"
#include "Material.h"
#include "DirectionalLight.h"
#include "PointLight.h"
#include "Light.h"
#include "LaplaceEigen.h"

struct ObjectPointerNode {
	Object* curr;
	ObjectPointerNode* next;
};

class Window
{
public:
	static ObjectPointerNode* head,* tail;
	static void addObjects(Object* toAdd);
	// Initializes your shader program(s) and uniform variable locations
	static bool initializeProgram();
	// Initialize your OBJ objects here
	static bool initializeObjects();
	// Make sure to clean up any data you dynamically allocated here
	static void cleanUp();
	// Creates a GLFW window context
	static GLFWwindow* createWindow(int width, int height);
	// Is called whenever the window is resized
	static void resizeCallback(GLFWwindow* window, int width, int height);
	// Is called on idle.
	static void idleCallback();
	// This renders to the glfw window. Add draw calls here
	static void displayCallback(GLFWwindow*);
	// Add your key press event handling here
	static void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods);
	static void cursorPosCallback(GLFWwindow* window, double xPos, double yPos);
	static void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods);
	static void scrollCallback(GLFWwindow* window, double xOffset, double yOffset);
	static glm::vec3 trackBallMapping(double xPos, double yPos);
};

#endif
