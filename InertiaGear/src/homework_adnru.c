#include "fem.h"


#ifndef NORHOSTEEL
double inertiaGearSteelRho()
{

//
// Modifier la valeur pour avoir la masse volumique de l'acier [kg/m3]
// Une tolerance de 10% sur la valeur est admise
//

    double rho = 7800;
    return rho;
}
#endif

#ifndef NOINERTIA
double inertiaGearInertia(femMesh *theMesh, femIntegration *theRule, double rho)
{
    double I = 0;
    double* X = theMesh->X;
    double* Y = theMesh->Y;
    int * elem = theMesh->elem;
    int nNode = theMesh->nLocalNode;

    double jacobian;
    double somme;
    
    for (size_t i = 0; i < theMesh->nElem; i++){
        jacobian = fabs(
            (X[elem[nNode*i + 1]] - X[elem[nNode*i + 0]]) * (Y[elem[nNode*i + 2]] - Y[elem[nNode*i + 0]])
            - (X[elem[nNode*i + 2]] - X[elem[nNode*i + 0]]) * (Y[elem[nNode*i + 1]] - Y[elem[nNode*i + 0]])
        );

        double xLoc_j;
        double yLoc_j;
        somme = 0;
        for (size_t j = 0; j < theRule->n; j++){
            xLoc_j = X[elem[nNode*i + 0]] * (1- theRule->xsi[j] - theRule->eta[j]) 
                + X[elem[nNode*i + 1]] * theRule->xsi[j] + X[elem[nNode*i + 2]] * theRule->eta[j];
            yLoc_j = Y[elem[nNode*i + 0]] * (1- theRule->xsi[j] - theRule->eta[j]) 
                + Y[elem[nNode*i + 1]] * theRule->xsi[j] + Y[elem[nNode*i + 2]] * theRule->eta[j];
            somme += theRule->weight[j] * ( xLoc_j*xLoc_j + yLoc_j*yLoc_j);
        }
        I += somme * jacobian * pow(10,-13);
    }
    
    return I * rho;
}
#endif

#ifndef NOREAD
femMesh *inertiaGearMeshRead(const char *filename)
{
    femMesh *theMesh = malloc(sizeof(femMesh));

    int i,trash;
    
    FILE* file = fopen(filename,"r");
    if (file == NULL) Error("No mesh file !");

    ErrorScan(fscanf(file, "Number of nodes %d \n", &theMesh->nNode));
    theMesh->X = malloc(sizeof(double)*theMesh->nNode);
    theMesh->Y = malloc(sizeof(double)*theMesh->nNode);
    
    for (i = 0; i < theMesh->nNode; ++i) {
        ErrorScan(fscanf(file,"%d : %le %le \n",&trash,&theMesh->X[i],&theMesh->Y[i])); }

    

    //ADDED
    theMesh->nLocalNode = 3;
    ErrorScan(fscanf(file, "Number of triangles %d \n", &theMesh->nElem));
    theMesh->elem = malloc(sizeof(int)*theMesh->nElem * theMesh->nLocalNode);
    //printf("nElem:%d\n",theMesh->nElem);
    for (i = 0; i < theMesh->nElem; ++i) {
        ErrorScan(fscanf(file,"%d : %d %d %d \n",
            &trash,&theMesh->elem[i*3],&theMesh->elem[i*3 + 1],&theMesh->elem[i*3 + 2])); 
        //printf("%d : %d %d %d \n",trash,theMesh->elem[i*3],theMesh->elem[i*3 + 1],theMesh->elem[i*3 + 2]);
    }
    //ADDED

    fclose(file);
    return theMesh;
}
#endif

#ifndef NOFREE
void inertiaGearMeshFree(femMesh *theMesh)
{
    free(theMesh->elem);
    free(theMesh->X);
    free(theMesh->Y);
    free(theMesh);
}
#endif

