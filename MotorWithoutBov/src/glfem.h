/*
 *  glfem.h
 *  Library for LEPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2017 UCL-EPL : Vincent Legat
 *  All rights reserved.
 *
 *  Pour GLFW (version utilisée 3.1)
 *  Pour l'installation de la librairie, voir http://www.glfw.org/
 *
 */

#ifndef _GLFEM_H_
#define _GLFEM_H_

#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>

#include "fem.h"

static float GLFEM_BLACK[4]         = {0.0,0.0,0.0,1.0};
static float GLFEM_BLUE[4]          = {0.0,0.0,1.0,1.0};
static float GLFEM_RED[4]           = {1.0,0.0,0.0,1.0};
static float GLFEM_GREEN[4]         = {0.0,1.0,0.0,1.0};
static float GLFEM_BACKGROUND[4]    = {0.9,0.9,0.8,0.0};



void        glfemDrawColorElement(float *x, float *y, double *u, int n);
void 		    glfemDrawElement(float *x, float *y, int n);
void 		    glfemDrawNodes(double* x, double* y,int n);
char        glfemGetAction();

void 		    glfemReshapeWindows(femMesh *theMesh, int width, int heigh);
void 		    glfemPlotSolution(femMesh *theMesh, double *u);
void 		    glfemPlotMesh(femMesh *theMesh);
void 		    glfemPlotEdges(femEdges *theEdges);
void 		    glfemPlotBnd(femEdges *theEdges);
void        glfemSetColor(float color[4]);
void        glfemSetScale(femMesh *theMesh, double *u);




void 		    glfemMessage(char *message);
void 		    glfemDrawMessage(int h, int v, char *message);
void 		    glfemSetRasterSize(int width, int height);
GLFWwindow* glfemInit(char *windowName);

static void glfemKeyCallback(GLFWwindow* self,int key,int scancode,int action,int mods);




#endif