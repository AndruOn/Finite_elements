
/*
 *  fem.h
 *  Library for LEPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2021 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#ifndef _FEM_H_
#define _FEM_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <assert.h>


#define ErrorScan(a)   femErrorScan(a,__LINE__,__FILE__)
#define Error(a)       femError(a,__LINE__,__FILE__)
#define Warning(a)     femWarning(a,  __LINE__, __FILE__)
#define FALSE 0 
#define TRUE  1



typedef enum {FEM_TRIANGLE,FEM_QUAD} femElementType; //!< Type d'élements du maillage (triangulaire ou quadrilataire)
typedef enum {FEM_FULL,FEM_BAND,FEM_ITER} femSolverType;
typedef enum {FEM_NO,FEM_XNUM,FEM_YNUM} femRenumType;

typedef struct {
    int *elem;
    double *X;
    double *Y;
    int nElem;
    int nNode;
    int nLocalNode;
    int *number;
} femMesh;

/**
 * Contient les paramètres de la méthode d'integration.
 * - n: Nombre de points
 * - xsi
 * - eta
 * - weight
 */
typedef struct {
    int n;                  //!< Nombre de points pour L'integration
    const double *xsi;      //!< Liste de taille n. Contenant les xsi
    const double *eta;      //!< Liste de taille n. Contenant les eta
    const double *weight;   //!< Liste de taille n. Contenant les weights
} femIntegration;

typedef struct {
    int elem[2];
    int node[2];
} femEdge;

typedef struct {
    femMesh *mesh;
    femEdge *edges;
    int nEdge;
    int nBoundary;
} femEdges;

typedef struct {
    int n;
    void (*x2)(double *xsi, double *eta);
    void (*phi2)(double xsi, double eta, double *phi);
    void (*dphi2dx)(double xsi, double eta, double *dphidxsi, double *dphideta);
} femDiscrete;

/**
 * Système complet:
 * - A : Matrice de raideur (sizexsize)
 * - B : Vecteur force (size)
 * - size
 */
typedef struct {
    double *B;
    double **A;
    int size;
} femFullSystem;

typedef struct {
    double *B;
    double **A;        
    int size;
    int band;        
} femBandSystem;

typedef struct {
    double *R;
    double *D;
    double *S;
    double *X; 
    double error;      
    int size;
    int iter;        
} femIterativeSolver;

/**
 * Contient:
 *  - le type de solveur
 *  - Le système complet (Matrice de raideur et vecteur force)
 */
typedef struct {
    femSolverType type;     //!< Type de solveur utilisé (BAND, FULL, ITER)
    femFullSystem *local;   //!< Contient une structure du système complet pour un élément local
    void *solver;           //!< Poite vers la fonction utiliser pour le solve le système
} femSolver;


/**
 * Contient toute les propiétés pour résoudre le systeme:
 * - mesh
 * - Règles d'espace
 * - Règles d'integration
 * - Solveur à utiliser
 * - Taille 
 * - Noeuds locaux
 * - Dirichlet
 * - Soluce
 * - Terme source
 * - Terme dircihlet
 */
typedef struct {
    femMesh *mesh;           //!< Mesh du probème 
    femDiscrete *space;      //!< Espace utilisé
    femIntegration *rule;    //!< Règle d'integration
    femSolver *solver;       //!< Solveur du système
    int size;                //!< Nombre de noeud du système
    int sizeLoc;             //!< Nombre de noeud par élément
    int *dirichlet;          //!< Tableau de taille size. Contient 0 ou 1 si le ième noeud est soumis ,ou non, à la condiiton de dircihlet. L'indice i correspond au ième noeud.
    double *soluce;          //!< Tableau de taille size. Contient la solution du problème.
    double sourceValue;      //!< Terme source
    double dirichletValue;   //!< Valeur condition de dirichlet
} femDiffusionProblem;


femIntegration      *femIntegrationCreate(int n, femElementType type);
void                 femIntegrationFree(femIntegration *theRule);

femMesh             *femMeshRead(const char *filename);
void                 femMeshWriteArray(const double *array, const int n, const char *filename);
void                 femMeshWrite(const femMesh* myMesh, const char *filename);
void                 femMeshFree(femMesh *theMesh);
void                 femMeshMakeCounterClockWise(femMesh *mesh);
void                 femMeshRenumber(femMesh *theMesh, femRenumType renumType);
int                  femMeshComputeBand(femMesh *theMesh);



femEdges*            femEdgesCreate(femMesh *theMesh);
void                 femEdgesFree(femEdges *theEdges);
void                 femEdgesPrint(femEdges *theEdges);
int                  femEdgesCompare(const void *edgeOne, const void *edgeTwo);

femDiscrete*         femDiscreteCreate(int n, femElementType type);
void                 femDiscreteFree(femDiscrete* mySpace);
void                 femDiscretePrint(femDiscrete* mySpace);
void                 femDiscreteXsi2(femDiscrete* mySpace, double *xsi, double *eta);
void                 femDiscretePhi2(femDiscrete* mySpace, double xsi, double eta, double *phi);
void                 femDiscreteDphi2(femDiscrete* mySpace, double xsi, double eta, double *dphidxsi, double *dphideta);

femSolver*           femSolverFullCreate(int size, int sizeLoc);
femSolver*           femSolverBandCreate(int size, int sizeLoc, int band);
femSolver*           femSolverIterativeCreate(int size, int sizeLoc);
void                 femSolverFree(femSolver* mySolver);
void                 femSolverInit(femSolver* mySolver);
void                 femSolverPrint(femSolver* mySolver);
void                 femSolverPrintInfos(femSolver* mySolver);
double*              femSolverEliminate(femSolver* mySolver);
void                 femSolverAssemble(femSolver* mySolver, double *Aloc, double *Bloc, double *Uloc, int *map, int nLoc);
double               femSolverGet(femSolver* mySolver, int i, int j);
int                  femSolverConverged(femSolver *mySolver);

femFullSystem*       femFullSystemCreate(int size);
void                 femFullSystemFree(femFullSystem* mySystem);
void                 femFullSystemInit(femFullSystem* mySystem);
void                 femFullSystemPrint(femFullSystem* mySystem);
void                 femFullSystemPrintInfos(femFullSystem* mySystem);
double*              femFullSystemEliminate(femFullSystem* mySystem);
void                 femFullSystemConstrain(femFullSystem* mySystem, int myNode, double value);
void                 femFullSystemAssemble(femFullSystem* mySystem, double *Aloc, double *Bloc, int *map, int nLoc);
double               femFullSystemGet(femFullSystem* mySystem, int i, int j);

femBandSystem*       femBandSystemCreate(int size, int band);
void                 femBandSystemFree(femBandSystem* myBandSystem);
void                 femBandSystemInit(femBandSystem *myBand);
void                 femBandSystemPrint(femBandSystem *myBand);
void                 femBandSystemPrintInfos(femBandSystem *myBand);
double*              femBandSystemEliminate(femBandSystem *myBand);
void                 femBandSystemAssemble(femBandSystem* myBandSystem, double *Aloc, double *Bloc, int *map, int nLoc);
double               femBandSystemGet(femBandSystem* myBandSystem, int i, int j);

femIterativeSolver*  femIterativeSolverCreate(int size);
void                 femIterativeSolverFree(femIterativeSolver* mySolver);
void                 femIterativeSolverInit(femIterativeSolver* mySolver);
void                 femIterativeSolverPrint(femIterativeSolver* mySolver);
void                 femIterativeSolverPrintInfos(femIterativeSolver* mySolver);
double*              femIterativeSolverEliminate(femIterativeSolver* mySolver);
void                 femIterativeSolverAssemble(femIterativeSolver* mySolver, double *Aloc, double *Bloc, double *Uloc, int *map, int nLoc);
double               femIterativeSolverGet(femIterativeSolver* mySolver, int i, int j);
int                  femIterativeSolverConverged(femIterativeSolver *mySolver);

femDiffusionProblem *femDiffusionCreate(const char *filename, femSolverType solverType, femRenumType renumType);
void                 femDiffusionFree(femDiffusionProblem *theProblem);
void                 femDiffusionMeshLocal(const femDiffusionProblem *theProblem, const int i, 
                              int *map, int *dirichlet, double *x, double *y, double *u);
void                 femDiffusionCompute(femDiffusionProblem *theProblem);
void                 femDiffusionRenumber(femDiffusionProblem *theProblem, femRenumType renumType);
int                  femDiffusionComputeBand(femDiffusionProblem *theProblem);

double               femMin(double *x, int n);
double               femMax(double *x, int n);
void                 femError(char *text, int line, char *file);
void                 femErrorScan(int test, int line, char *file);
void                 femWarning(char *text, int line, char *file);
#endif
