#include "Window.h"

std::string bunnyFilename = ".\\Objects\\bunny.obj";
std::string sphereFilename = ".\\Objects\\sphere.obj";
ObjectPointerNode* Window::head, * Window::tail;
PointCloud* bunny;
Light* directionalLight, * pointLight;
/*
 * Declare your variables below. Unnamed namespace is used here to avoid
 * declaring global or static variables.
 */
namespace
{
	int width, height;
	std::string windowTitle("GLFW Starter Project");

	Cube* cube;
	PointCloud* objectPoints;
	Object* currentObj, *controlledObject;
	Light* currentLight;

	glm::vec3 eye(0, 0, 20); // Camera position.
	glm::vec3 center(0, 0, 0); // The point we are looking at.
	glm::vec3 up(0, 1, 0); // The up direction of the camera.
	float fovy = 60;
	float near = 1;
	float far = 1000;
	glm::mat4 view = glm::lookAt(eye, center, up); // View matrix, defined by eye, center and up.
	glm::mat4 projection; // Projection matrix.

	GLuint currentProgram, normalProgram, lightProgram; // The shader program id.
	const int PROGRAM_COUNT = 2;
	GLuint programs[PROGRAM_COUNT];
	int programInd = 0;

	GLuint projectionLoc; // Location of projection in shader.
	GLuint viewLoc; // Location of view in shader.
	GLuint modelLoc; // Location of model in shader.
	GLuint colorLoc; // Location of color in shader.
	GLuint projectionTILoc, viewTILoc, modelTILoc; // Location of transpose inverse
	GLuint positionLoc; //Location of light source
	GLuint eyeLoc; // Location of eye in the world
	GLuint dirLightOnLoc;
	GLuint pointLightOnLoc;

	bool mouseLeftPressed = false;
	bool mouseRightPressed = false;
	double prevX, prevY;
	bool dirLightOn, pointLightOn;

	// config for LE
	int PARTICLE_NUM = 100;
	int GRID_RES = 100;
	int DENS_RES = 100;
	int BASIS_DIM = 16;
	std::vector<PointCloud*> particles;
	LaplaceEigen* le_object;
};

void setShaderLoc() {
	glUseProgram(currentProgram);
	// Get the locations of uniform variables.
	projectionLoc   = glGetUniformLocation(currentProgram, "projection");
	viewLoc         = glGetUniformLocation(currentProgram, "view");
	modelLoc        = glGetUniformLocation(currentProgram, "model");
	colorLoc        = glGetUniformLocation(currentProgram, "color");
	projectionTILoc = glGetUniformLocation(currentProgram, "projectionTI");
	viewTILoc       = glGetUniformLocation(currentProgram, "viewTI");
	modelTILoc      = glGetUniformLocation(currentProgram, "modelTI");
	positionLoc     = glGetUniformLocation(currentProgram, "lightPosition");
	eyeLoc          = glGetUniformLocation(currentProgram, "eye");
	dirLightOnLoc   = glGetUniformLocation(currentProgram, "dirLightOn");
	pointLightOnLoc = glGetUniformLocation(currentProgram, "pointLightOn");
}

bool Window::initializeProgram()
{
	// Create a shader program with a vertex shader and a fragment shader.
	normalProgram = LoadShaders("shaders/shader.vert", "shaders/shader.frag");
	lightProgram = LoadShaders("shaders/lightShader.vert", "shaders/lightShader.frag");
	programs[0] = normalProgram;
	programs[1] = lightProgram;
	
	// Check the shader programs.
	if (!normalProgram || !lightProgram)
	{
		std::cerr << "Failed to initialize shader program" << std::endl;
		return false;
	}

	currentProgram = normalProgram;
	setShaderLoc();
	return true;
}

void Window::addObjects(Object* toAdd) {
	if (head == nullptr) {
		head = new ObjectPointerNode();
		head->curr = toAdd;
		tail = head;
	}

	tail->next = new ObjectPointerNode();
	tail->curr = toAdd;
	tail = tail->next;
}

bool Window::initializeObjects()
{
	// Create a point cloud consisting of bunny
	objectPoints = new PointCloud(bunnyFilename, 1);
	objectPoints->setMaterial(Material::DefinedMaterial::REGULAR, currentProgram);
	bunny = objectPoints;
	addObjects(objectPoints);

	// Create light
	directionalLight = new DirectionalLight(currentProgram, 3.f, 3.f, 3.f);
	pointLight = new PointLight(.0f, 1.f, 1.f);

	controlledObject = currentObj = bunny;
	currentLight = directionalLight;

	le_object = new LaplaceEigen(GRID_RES, BASIS_DIM, DENS_RES);
	le_object->add_particles(PARTICLE_NUM);

	for (int i = 0; i < BASIS_DIM; i++) {
		le_object->forces_dw[i] = 10 * ((double)std::rand()) / RAND_MAX;
	}

	for (int i = 0; i < le_object->particles.size(); i++) {
		PointCloud* particle = new PointCloud(sphereFilename, 1, 5*(2 * le_object->particles[i][0]-1), 5*(2*le_object->particles[i][1]-1));
		particle->scale(0.3, 0.3, 0.3);
		particle->translate();
		particles.push_back(particle);
	}

	return true;
}

void Window::cleanUp()
{
	// Deallcoate the objects.
	delete cube;
	delete objectPoints;

	// Delete the shader program.
	glDeleteProgram(normalProgram);
	// glDeleteProgram(lightProgram);
}

GLFWwindow* Window::createWindow(int width, int height)
{
	// Initialize GLFW.
	if (!glfwInit())
	{
		std::cerr << "Failed to initialize GLFW" << std::endl;
		return NULL;
	}

	// 4x antialiasing.
	glfwWindowHint(GLFW_SAMPLES, 4);

#ifdef __APPLE__ 
	// Apple implements its own version of OpenGL and requires special treatments
	// to make it uses modern OpenGL.

	// Ensure that minimum OpenGL version is 3.3
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	// Enable forward compatibility and allow a modern OpenGL context
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

	// Create the GLFW window.
	GLFWwindow* window = glfwCreateWindow(width, height, windowTitle.c_str(), NULL, NULL);

	// Check if the window could not be created.
	if (!window)
	{
		std::cerr << "Failed to open GLFW window." << std::endl;
		glfwTerminate();
		return NULL;
	}

	// Make the context of the window.
	glfwMakeContextCurrent(window);

#ifndef __APPLE__
	// On Windows and Linux, we need GLEW to provide modern OpenGL functionality.

	// Initialize GLEW.
	if (glewInit())
	{
		std::cerr << "Failed to initialize GLEW" << std::endl;
		return NULL;
	}
#endif

	// Set swap interval to 1.
	glfwSwapInterval(0);


	// Call the resize callback to make sure things get drawn immediately.
	Window::resizeCallback(window, width, height);

	return window;
}

void Window::resizeCallback(GLFWwindow* window, int w, int h)
{
#ifdef __APPLE__
	// In case your Mac has a retina display.
	glfwGetFramebufferSize(window, &width, &height);
#endif
	width = w;
	height = h;

	// Set the viewport size.
	glViewport(0, 0, width, height);

	// Set the projection matrix.
	projection = glm::perspective(glm::radians(fovy),
		(float)width / (float)height, near, far);
}

void Window::idleCallback()
{
	// Perform any updates as necessary.
	ObjectPointerNode* current = head;
	while (current->curr) {
		current->curr->update();
		current = current->next;
	}
}

void Window::displayCallback(GLFWwindow* window)
{
	// Switch back to using OpenGL's rasterizer
	glUseProgram(currentProgram);
	// Clear the color and depth buffers.
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Specify the values of the uniform variables we are going to use.
	glm::mat4 model = currentObj->getModel();
	glm::vec3 color = currentObj->getColor();
	glUniformMatrix4fv(projectionLoc, 1, GL_FALSE, glm::value_ptr(projection));
	glUniformMatrix4fv(viewLoc, 1, GL_FALSE, glm::value_ptr(view));
	glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));

	glm::mat4 modelTI = glm::inverse(glm::transpose(model));
	glm::mat4 projectionTI = glm::inverse(glm::transpose(projection));
	glm::mat4 viewTI = glm::inverse(glm::transpose(view));
	glUniformMatrix4fv(projectionTILoc, 1, GL_FALSE, glm::value_ptr(projectionTI));
	glUniformMatrix4fv(viewTILoc, 1, GL_FALSE, glm::value_ptr(viewTI));
	glUniformMatrix4fv(modelTILoc, 1, GL_FALSE, glm::value_ptr(modelTI));

	glUniform3fv(eyeLoc, 1, glm::value_ptr(eye));
	glUniform3fv(colorLoc, 1, glm::value_ptr(color));

	// Send information about material to binded shader
	((PointCloud*)currentObj)->getMaterial()->setShaderProgram(currentProgram);
	((PointCloud*)currentObj)->getMaterial()->sendToShader();

	if (dirLightOn) {
		directionalLight->setShaderProgram(currentProgram);
		directionalLight->sendToShader();
		directionalLight->update();
	}

	if (pointLightOn) {
		pointLight->setShaderProgram(currentProgram);
		pointLight->sendToShader();
	}
	// Render the object.
	currentObj->draw();

	le_object->step();
	le_object->advect_particles();
	for (int i = 0; i < PARTICLE_NUM; i++) {
		delete particles[i];
		particles[i] = new PointCloud(sphereFilename, 1, 5 * (2 * le_object->particles[i][0] - 1), 5 * (2 * le_object->particles[i][1] - 1));
		particles[i]->scale(0.3, 0.3, 0.3);
		glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(particles[i]->getModel()));
		particles[i]->draw();
	}

	if (pointLightOn) {
		// This should not be drawn before other objects because it resets model for itself
		((PointLight*)pointLight)->showRepresentation(modelLoc, colorLoc);
	}

	// Gets events, including input such as keyboard and mouse or window resizing.
	glfwPollEvents();
	// Swap buffers.
	glfwSwapBuffers(window);

}

void Window::keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	// Check for a key press.
	ObjectPointerNode* current;
	if (action == GLFW_PRESS)
	{
		switch (key)
		{
		case GLFW_KEY_ESCAPE:
			// Close the window. This causes the program to also terminate.
			glfwSetWindowShouldClose(window, GL_TRUE);
			break;

		case 290:
			currentObj = bunny;
			controlledObject = currentObj;
			break;

		case GLFW_KEY_0:
			controlledObject = (controlledObject == currentObj) ? 
				((PointLight*)pointLight)->representation : currentObj;
			break;

		case GLFW_KEY_1:
			if (currentProgram != lightProgram) break;
			glUniform1i(dirLightOnLoc, (GLuint)1);
			glUniform1i(pointLightOnLoc, (GLuint)0);
			dirLightOn = true;
			pointLightOn = false;
			controlledObject = currentObj;
			break;

		case GLFW_KEY_2:
			if (currentProgram != lightProgram) break;
			glUniform1i(dirLightOnLoc, (GLuint)1);
			glUniform1i(pointLightOnLoc, (GLuint)1);
			dirLightOn = true;
			pointLightOn = true;
			break;

		case GLFW_KEY_3:
			if (currentProgram != lightProgram) break;
			break;

		case GLFW_KEY_4:
			if (currentProgram != lightProgram) break;
			dirLightOn = ! dirLightOn;
			glUniform1i(dirLightOnLoc, (GLuint) dirLightOn);
			break;

		case 'N':
			programInd = (programInd + 1) % PROGRAM_COUNT;
			currentProgram = programs[programInd];

			setShaderLoc();
			glUniform1i(dirLightOnLoc, (GLuint)1);
			glUniform1i(pointLightOnLoc, (GLuint)0);
			dirLightOn = true;
			pointLightOn = false;
			break;

		case 'P':
			current = head;
			while (current->curr) {
				PointCloud* pcObject = (PointCloud*)current->curr;
				if (mods) {
					pcObject->updatePointSize(pcObject->getPointSize() + 1);
				}
				else {
					pcObject->updatePointSize(pcObject->getPointSize() - 1);
				}
				current = current->next;
			}
			break;

		case 'A':
			current = head;
			while (current->curr) {
				PointCloud* pcObject = (PointCloud*)current->curr;
				pcObject->translate(-1, 0, 0);
				current = current->next;
			}
			break;

		case 'D':
			current = head;
			while (current->curr) {
				PointCloud* pcObject = (PointCloud*)current->curr;
				pcObject->translate(1, 0, 0);
				current = current->next;
			}
			break;

		case 'W':
			current = head;
			while (current->curr) {
				PointCloud* pcObject = (PointCloud*)current->curr;
				pcObject->translate(0, 1, 0);
				current = current->next;
			}
			break;

		case 'X':
			current = head;
			while (current->curr) {
				PointCloud* pcObject = (PointCloud*)current->curr;
				pcObject->translate(0, -1, 0);
				current = current->next;
			}
			break;

		case 'Z':
			current = head;
			while (current->curr) {
				PointCloud* pcObject = (PointCloud*)current->curr;
				if (mods) {
					pcObject->translate(0, 0, 1);
				}
				else {
					pcObject->translate(0, 0, -1);
				}
				current = current->next;
			}
			break;

		case 'S':
			current = head;
			while (current->curr) {
				PointCloud* pcObject = (PointCloud*)current->curr;
				if (mods) {
					pcObject->scale(1.1, 1.1, 1.1);
				}
				else {
					pcObject->scale(0.9, 0.9, 0.9);
				}
				current = current->next;
			}
			break;

		case 'R':
			current = head;
			while (current->curr) {
				PointCloud* pcObject = (PointCloud*)current->curr;
				if (mods) {
					pcObject->cancelScaleAndRot();
				}
				else {
					pcObject->cancelTranslate();
				}
				current = current->next;
			}
			break;

		default:
			break;
		}
	}
}

void Window::cursorPosCallback(GLFWwindow* window, double xPos, double yPos)
{
	// Based on movement of cursor positions, generate rotate axis or translation axis
	if (mouseLeftPressed) {
		glm::vec3 previousPos = trackBallMapping(prevX, prevY);
		glm::vec3 currentPos = trackBallMapping(xPos, yPos);
		glm::vec3 rotateAxis = glm::vec4(glm::cross(previousPos, currentPos), .0f);
		prevX = xPos;
		prevY = yPos;

		float radianAngle = glm::acos( std::min(std::abs( glm::dot(previousPos, currentPos) ), 1.0f) );
		((PointCloud*) controlledObject)->rotate(radianAngle, rotateAxis);
	}

	if (mouseRightPressed) {
		// Translate x & y
		((PointCloud*) controlledObject)->translate((xPos - prevX) * 25 / width, (prevY - yPos) * 25 / height, 0);
		prevX = xPos;
		prevY = yPos;
	}
}

void Window::mouseButtonCallback(GLFWwindow* window, int button, int action, int mods)
{
	// Set mouseLeftPressed / mouseRightPressed when being pressed
	// Record prevX and prevY for mouseLeftPressed with initial positions
	if (button == GLFW_MOUSE_BUTTON_LEFT){
		if (action == GLFW_PRESS) {
			mouseLeftPressed = true;

			glfwGetCursorPos(window, &prevX, &prevY);
		}
		else if (action == GLFW_RELEASE) {
			mouseLeftPressed = false;
		}
	}
	else if (button == GLFW_MOUSE_BUTTON_RIGHT) {
		if (action == GLFW_PRESS) {
			mouseRightPressed = true;

			glfwGetCursorPos(window, &prevX, &prevY);
		}
		else if (action == GLFW_RELEASE) {
			mouseRightPressed = false;
		}
	}
}

void Window::scrollCallback(GLFWwindow* window, double xOffset, double yOffset)
{
	if (yOffset > 0) {
		((PointCloud*)currentObj)->translate(.0f, .0f, .3f);
	}
	else if (yOffset < 0) {
		((PointCloud*)currentObj)->translate(.0f, .0f, -.3f);
	}
}

glm::vec3 Window::trackBallMapping(double xPos, double yPos)
{
	// Need a normalized result
	glm::vec3 result(.0f, .0f, .0f);
	float d;
	result[0] = (2.0 * xPos - width) / width;
	result[1] = (height - 2.0 * yPos) / height;
	d = result.length();
	d = (d < 1.0f) ? d : 1.0f;
	result[2] = sqrtf(1.001 - ((double)d) * d);
	
	return glm::normalize(result);
}