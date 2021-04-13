
#include"fem.h"


double convergenceSource(double x, double y)
{

//
// A obtenir avec sympy :-)
//
    double f = 1.0;
    return f; 
}

double convergenceSoluce(double x, double y, double *u)
{
    double a = sqrt(2.0);
    u[0] = x*y*(1-x)*(1-y)*atan(20*(x+y)/a - 16);
    u[1] = y*(y - 1)*(10*a*x*(x - 1)
      + (2*x - 1)*(4*pow(5*a*(x + y) - 8, 2) + 1)*atan(10*a*(x + y) - 16))
      /(4*pow(5*a*(x + y) - 8, 2) + 1);
    u[2] = x*(x - 1)*(10*a*y*(y - 1) 
      + (2*y - 1)*(4*pow(5*a*(x + y) - 8, 2) + 1)*atan(10*a*(x + y) - 16))
      /(4*pow(5*a*(x + y) - 8, 2) + 1);
    return u[0];
}

double convergenceEstimateRate(double *errors, int n, double ratio)
{

     
     double rate = 0;
     
//
//     A modifier
//     
     
     return rate;
}


void femDiffusionComputeError(femDiffusionProblem *theProblem, 
                                    double(*soluce)(double,double,double*))
{


//  
//     femMesh *theMesh = theProblem->mesh;
//     femIntegration *theRule = theProblem->rule;
//     femDiscrete *theSpace = theProblem->space;
// 
// 
//      A modifier
//
//       
    
    theProblem->errorSoluceL2 = 0.0;
    theProblem->errorSoluceH1 = 0.0;
    theProblem->errorInterpolationL2 = 0.0;
    theProblem->errorInterpolationH1 = 0.0;
    
  
}




femMesh *femMeshCreateBasicSquare(int n) 
{
    femMesh *theMesh = malloc(sizeof(femMesh));

    // Toujours le maillage élémentaire :-)
    // A modifier [begin]

    n = 2;   
        
    theMesh->nNode = (n+1)*(n+1);
    theMesh->X = malloc(sizeof(double)*theMesh->nNode);
    theMesh->Y = malloc(sizeof(double)*theMesh->nNode);
    double *X = theMesh->X;
    double *Y = theMesh->Y;
    X[0]= 0.0; Y[0]= 0.0;
    X[1]= 0.0; Y[1]= 0.5;
    X[2]= 0.0; Y[2]= 1.0;
    X[3]= 0.5; Y[3]= 0.0;
    X[4]= 0.5; Y[4]= 0.5;
    X[5]= 0.5; Y[5]= 1.0;
    X[6]= 1.0; Y[6]= 0.0;
    X[7]= 1.0; Y[7]= 0.5;
    X[8]= 1.0; Y[8]= 1.0;
    
    theMesh->nElem = 2*n*n;
    theMesh->nLocalNode = 3;
    theMesh->elem = malloc(sizeof(int)*3*theMesh->nElem);
    int *elem = theMesh->elem;
    int i;
    i=0; elem[3*i] = 0; elem[3*i+1] = 4; elem[3*i+2] = 1;
    i++; elem[3*i] = 3; elem[3*i+1] = 4; elem[3*i+2] = 0;
    i++; elem[3*i] = 1; elem[3*i+1] = 5; elem[3*i+2] = 2;
    i++; elem[3*i] = 4; elem[3*i+1] = 5; elem[3*i+2] = 1;
    i++; elem[3*i] = 3; elem[3*i+1] = 7; elem[3*i+2] = 4;
    i++; elem[3*i] = 6; elem[3*i+1] = 7; elem[3*i+2] = 3;
    i++; elem[3*i] = 4; elem[3*i+1] = 8; elem[3*i+2] = 5;
    i++; elem[3*i] = 7; elem[3*i+1] = 8; elem[3*i+2] = 4;
    
    // A modifier [end]

    theMesh->number = malloc(sizeof(int)*theMesh->nNode); 
    for (int i = 0; i < theMesh->nNode; i++) 
          theMesh->number[i] = i;     
    return theMesh;
}




