#include "fem.h"


#ifndef NORHOSTEEL
double inertiaGearSteelRho()
{

//
// Modifier la valeur pour avoir la masse volumique de l'acier [kg/m3]
// Une tolerance de 10% sur la valeur est admise
//

    double rho = 7800.0;
    return rho;
}
#endif


//typedef struct {
  //  int n;
    //const double *xsi;
    //const double *eta;
    //const double *weight;
//} femIntegration;


#ifndef NOINERTIA
double inertiaGearInertia(femMesh *theMesh, femIntegration *theRule, double rho)
{
    double I = 0.0;

    
    //double J = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);

    double Jacobien; 
    double posX[3]; 
    double posY[3]; 
    double sum = 0;
    
    
    for (int i = 0; i < theMesh->nElem; i++)
    {
        
        int T[3] = {theMesh -> elem[3*i], theMesh -> elem[3*i+1], theMesh -> elem[3*i+2]};
        posX[0] = theMesh->X[T[0]]; 
        posX[1] = theMesh->X[T[1]]; 
        posX[2] = theMesh->X[T[2]]; 
        posY[0] = theMesh->Y[T[0]]; 
        posY[1] = theMesh->Y[T[1]]; 
        posY[2] = theMesh->Y[T[2]]; 
        Jacobien = fabs((posX[1] - posX[0])*(posY[2]-posY[0]) - (posX[2] - posX[0])*(posY[1]-posY[0])); 

        for (size_t j = 0; j < theRule->n; j++)
        {
            double xLoc = posX[0]*(1-theRule->xsi[j]-theRule->eta[j]) + posX[1]*(theRule->xsi[j]) + posX[2]*theRule->eta[j];
            double yLoc = posY[0]*(1-theRule->xsi[j]-theRule->eta[j]) + posY[1]*(theRule->xsi[j]) + posY[2]*theRule->eta[j];
            sum += theRule->weight[j]*(xLoc*xLoc + yLoc*yLoc)*Jacobien;
        }
        


    }

    return rho*sum * pow(10,-13);
}
#endif

#ifndef NOREAD
femMesh *inertiaGearMeshRead(const char *filename)
{
    femMesh *theMesh = malloc(sizeof(femMesh));

    int i,trash,*elem;
    
    FILE* file = fopen(filename,"r");
    if (file == NULL) Error("No mesh file !");

    ErrorScan(fscanf(file, "Number of nodes %d \n", &theMesh->nNode));
    theMesh->X = malloc(sizeof(double)*theMesh->nNode);
    theMesh->Y = malloc(sizeof(double)*theMesh->nNode);
    for (i = 0; i < theMesh->nNode; ++i) {
        ErrorScan(fscanf(file,"%d : %le %le \n",&trash,&theMesh->X[i],&theMesh->Y[i])); }
    
    ErrorScan(fscanf(file, "Number of triangles %d \n", &theMesh->nElem));
    theMesh-> elem = malloc(sizeof(double)*theMesh->nElem*3);
    for (i = 0; i < theMesh->nElem; ++i) {
        ErrorScan(fscanf(file,"%d : %d %d %d \n",&trash,&theMesh->elem[3*i],&theMesh->elem[3*i+1], &theMesh->elem[3*i+2])); }




//
// A completer :-)
//     Allocation dynamique du tableau d'appartenance
//     Lecture du tableau d'appartenance
// 
// Le premiere partie fournie de la fonction peut largement servir d'inspiration....
// N'oubliez pas : Google est votre ami : par exemple google C malloc peut etre tres utile !
//     
//
    theMesh->nLocalNode = 3;
    
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

