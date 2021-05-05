#include "motor.h"



//je passe d'un motorMesh a un femmesh pour ne pas faire de cast deg dans femMeshRead
femMesh* motormesh_to_femmesh(motorMesh* theMotorMesh){
    
    femMesh* myFemMesh = malloc(sizeof(femMesh)); 
    myFemMesh->elem = theMotorMesh->elem; 
    myFemMesh->X = theMotorMesh->X; 
    myFemMesh->Y = theMotorMesh->Y; 
    myFemMesh->nElem = theMotorMesh->nElem; 
    myFemMesh->nNode = theMotorMesh->nNode; 
    myFemMesh->nLocalNode = theMotorMesh->nLocalNode; 
    myFemMesh->number = malloc(sizeof(int)* myFemMesh->nNode); 
    for (int i = 0; i < myFemMesh->nNode; i++) {myFemMesh->number[i] = i; }
    return myFemMesh; 
}


/*fonctions modifi�s de fem.c*/



void femDiffusionComputeConrad(femDiffusionProblem* theProblem, motor *theMotor)
{
    femMesh* theMesh = theProblem->mesh;
    femIntegration* theRule = theProblem->rule;
    femDiscrete* theSpace = theProblem->space;
    femSolver* theSolver = theProblem->solver;
    int* number = theProblem->mesh->number;

    double dirichlet = theProblem->dirichletValue;

    if (theSpace->n > 4) Error("Unexpected discrete space size !");
    double Xloc[4], Yloc[4], phi[4], dphidxsi[4], dphideta[4], dphidx[4], dphidy[4];
    double Uloc[4];
    int iEdge, iElem, iInteg, i, j, map[4], ctr[4];
    double** A = theSolver->local->A;
    double* Aloc = theSolver->local->A[0];
    double* Bloc = theSolver->local->B;

    motorMesh *themotorMesh = theMotor->mesh; 
   

    for (iElem = 0; iElem < theMesh->nElem; iElem++) {

        //je r�cup�re le domaine du moteur dans lequel on est occup�
        int domain = themotorMesh->domain[iElem]; 
        //je prends les constantes sp�cifiques au domaine 
        double js = theMotor->js[domain]; 
        double mu = theMotor->mu[domain];

        for (i = 0; i < theSpace->n; i++)      Bloc[i] = 0;
        for (i = 0; i < (theSpace->n) * (theSpace->n); i++) Aloc[i] = 0;
        femDiffusionMeshLocal(theProblem, iElem, map, ctr, Xloc, Yloc, Uloc);
        for (iInteg = 0; iInteg < theRule->n; iInteg++) {
            double xsi = theRule->xsi[iInteg];
            double eta = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];
            femDiscretePhi2(theSpace, xsi, eta, phi);
            femDiscreteDphi2(theSpace, xsi, eta, dphidxsi, dphideta);
            double dxdxsi = 0;
            double dxdeta = 0;
            double dydxsi = 0;
            double dydeta = 0;
            for (i = 0; i < theSpace->n; i++) {
                dxdxsi += Xloc[i] * dphidxsi[i];
                dxdeta += Xloc[i] * dphideta[i];
                dydxsi += Yloc[i] * dphidxsi[i];
                dydeta += Yloc[i] * dphideta[i];
            }
            double jac = dxdxsi * dydeta - dxdeta * dydxsi;
            for (i = 0; i < theSpace->n; i++) {
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac;
            }
            for (i = 0; i < theSpace->n; i++) {
                for (j = 0; j < theSpace->n; j++) {
                    //on multiplie chaque entr�e par 1/mu 
                    A[i][j] += (1.0/mu) * ((dphidx[i] * dphidx[j]
                        + dphidy[i] * dphidy[j]) * jac * weight);
                }
            }
            for (i = 0; i < theSpace->n; i++) {
                //j'ai remplac� le temre source par le terme source du domaine dans lequel on est 
                Bloc[i] += phi[i] * jac * js * weight;
            }
        }
        for (i = 0; i < theSpace->n; i++)
            if (ctr[i] == 1) femFullSystemConstrain(theSolver->local, i, dirichlet);
        femSolverAssemble(theSolver, Aloc, Bloc, Uloc, map, theSpace->n);
    }

    double* soluce = femSolverEliminate(theSolver);
    for (i = 0; i < theProblem->size; i++)
        theProblem->soluce[i] += soluce[number[i]];
}

//j'ai remplac� le file par theMotor (c'est l� qu'on a toute l'info)
//fflush(stdout) 
femDiffusionProblem* femDiffusionCreateConrad(motor* theMotor, femSolverType solverType, femRenumType renumType)
{
    int i, band;

    femDiffusionProblem* theProblem = malloc(sizeof(femDiffusionProblem));
    printf("ligne 85"); fflush(stdout);

    //plus besoin de cr�er le mesh � partir d'un fichier, il nous est donn� => je le r�cup�re
    //Je le change en femMesh au lieu de motMesh car les fonctions sont faites pour marcher sur des femMesh
    femMesh* myFemMesh  = motormesh_to_femmesh(theMotor->mesh); 
    printf("ligne 86"); fflush(stdout);

    theProblem->mesh = myFemMesh;
    printf("ligne 88"); fflush(stdout); 

    //plus besoin de prendre en compte le cas de nLocalnode = 4, on sait qu'on travaille avec des triangles
    theProblem->space = femDiscreteCreate(3, FEM_TRIANGLE);
    printf("ligne 99"); fflush(stdout);
    theProblem->rule = femIntegrationCreate(3, FEM_TRIANGLE);
    printf("ligne 102"); fflush(stdout);

    theProblem->size = theProblem->mesh->nNode;
    theProblem->sizeLoc = theProblem->mesh->nLocalNode;
    femMeshRenumber(theProblem->mesh, renumType);
    theProblem->sourceValue = 1.0;

    //je mets la fronti�re � 0 (conditions fronti�re)
    theProblem->dirichletValue = 0.0;
    printf("ligne 117"); fflush(stdout);


    theProblem->dirichlet = malloc(sizeof(int) * theProblem->size);
    for (i = 0; i < theProblem->size; i++)
        theProblem->dirichlet[i] = 0;
    femEdges* theEdges = femEdgesCreate(theProblem->mesh);
    for (i = 0; i < theEdges->nEdge; i++) {
        if (theEdges->edges[i].elem[1] < 0) {
            theProblem->dirichlet[theEdges->edges[i].node[0]] = 1;
            theProblem->dirichlet[theEdges->edges[i].node[1]] = 1;
        }
    }
    femEdgesFree(theEdges);
    printf("ligne 131"); fflush(stdout);
    switch (solverType) {
    case FEM_FULL:
        theProblem->solver = femSolverFullCreate(theProblem->size,
            theProblem->sizeLoc); break;
    case FEM_BAND:
        band = femMeshComputeBand(theProblem->mesh);
        theProblem->solver = femSolverBandCreate(theProblem->size,
            theProblem->sizeLoc, band); break;
    case FEM_ITER:
        theProblem->solver = femSolverIterativeCreate(theProblem->size,
            theProblem->sizeLoc); break;
    default: Error("Unexpected solver option");
    }

    theProblem->soluce = malloc(sizeof(double) * theProblem->size);
    for (i = 0; i < theProblem->size; i++)
        theProblem->soluce[i] = 0;
    printf("ligne 149"); fflush(stdout);


    return theProblem;
}


//j'enlève le free mesh
void femDiffusionFreeConrad(femDiffusionProblem *theProblem)
{
    femIntegrationFree(theProblem->rule);
    femDiscreteFree(theProblem->space);
    //femMeshFree(theProblem->mesh);
    femSolverFree(theProblem->solver);
    free(theProblem->dirichlet);
    free(theProblem->soluce);
    free(theProblem);
}


//
// ========= Projet � r�aliser ===================
//

void motorAdaptMesh(motor *theMotor, double delta)
{
    motorMesh *theMesh = theMotor->mesh;
    
    double x,y;
    for(int i = 0; i < theMesh->nNode; ++i){
        if  (theMotor->movingNodes[i] == 1){
            x = theMesh->X[i]*cos(delta) - theMesh->Y[i]*sin(delta);
            y = theMesh->X[i]*sin(delta) + theMesh->Y[i]*cos(delta);
            theMesh->X[i] = x;
            theMesh->Y[i] = y; }}
    theMotor->theta += delta;
}


double motorComputeCouple(motor *theMotor)
{
    return 0.0;

}

void motorComputeCurrent(motor *theMotor)
{
    return;
}

void motorComputeMagneticPotential(motor *theMotor)
{
    //je cr�� le probl�me
    femDiffusionProblem* myProblem = femDiffusionCreateConrad(theMotor, FEM_BAND, FEM_NO);
    printf("ligne 189"); fflush(stdout);
    //je solve le probl�me
    femDiffusionComputeConrad(myProblem, theMotor); 
    //je rentre les valeurs nodales du potentiel magn�tique dans la matrice correspondante
    theMotor->a = myProblem->soluce; 
    femDiffusionFreeConrad(myProblem); 

    return;
    
}

//
// ========= Projet � r�aliser ===================
//







