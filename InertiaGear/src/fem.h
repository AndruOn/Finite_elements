
/*
 *  fem.h
 *  Library for MECA1120 : Finite Elements for dummies
 *
 *  Copyright (C) 2016 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#ifndef _FEM_H_
#define _FEM_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define ErrorScan(a)   femErrorScan(a,__LINE__,__FILE__)
#define Error(a)       femError(a,__LINE__,__FILE__)
#define Warning(a)     femWarning(a,  __LINE__, __FILE__)
#define FALSE 0 
#define TRUE  1


typedef enum {FEM_TRIANGLE,FEM_QUAD} femElementType;

typedef struct {
    int *elem;
    double *X;
    double *Y;
    int nElem;
    int nNode;
    int nLocalNode;
} femMesh;


typedef struct {
    int n;
    const double *xsi;
    const double *eta;
    const double *weight;
} femIntegration;


femIntegration  *femIntegrationCreate(int n, femElementType type);
void             femIntegrationFree(femIntegration *theRule);

femMesh         *inertiaGearMeshRead(const char *filename);
double           inertiaGearSteelRho();
double           inertiaGearInertia(femMesh *theMesh, femIntegration *theRule, double rho);
void             inertiaGearMeshFree(femMesh *theMesh);



double           femMin(double *x, int n);
double           femMax(double *x, int n);
void             femError(char *text, int line, char *file);
void             femErrorScan(int test, int line, char *file);
void             femWarning(char *text, int line, char *file);

#endif
