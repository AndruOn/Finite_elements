#include "motor.h"
#include "fem.h"


motorMesh* theMesh;
motor* theMotor_;


//Q1
femMesh* motorMesh_to_femMesh(motorMesh *theMotorMesh){
    femMesh *thefemMesh = malloc(sizeof(femMesh));
    thefemMesh->elem = theMotorMesh->elem;
    thefemMesh->X = theMotorMesh->X;
    thefemMesh->Y = theMotorMesh->Y;
    thefemMesh->nElem = theMotorMesh->nElem;
    thefemMesh->nNode = theMotorMesh->nNode;
    thefemMesh->nLocalNode= theMotorMesh->nLocalNode;
    thefemMesh->number = malloc(sizeof(int)*thefemMesh->nNode); 
    for (int i = 0; i < thefemMesh->nNode; i++) 
          thefemMesh->number[i] = i; 
    return thefemMesh;
}

femDiffusionProblem *femDiffusionCreate_mesh(femMesh *thefemMesh, femSolverType solverType, femRenumType renumType)
{
    int i,band;
    
    femDiffusionProblem *theProblem = malloc(sizeof(femDiffusionProblem));
    theProblem->mesh  = thefemMesh;    //seule ligne changée       
    
    //tjrs avec des triangles
    theProblem->space = femDiscreteCreate(3,FEM_TRIANGLE);
    theProblem->rule = femIntegrationCreate(3,FEM_TRIANGLE);


    theProblem->size = theProblem->mesh->nNode;
    theProblem->sizeLoc = theProblem->mesh->nLocalNode;
    femMeshRenumber(theProblem->mesh,renumType);
    theProblem->sourceValue = 1.0;
    theProblem->dirichletValue = 0.0; // zero sur les bord
    
    
    theProblem->dirichlet = malloc(sizeof(int)*theProblem->size);
    for (i = 0; i < theProblem->size; i++) 
        theProblem->dirichlet[i] = 0;
    
    //faire une fois en global-----------------------
    //if(edges != NULL){
        femEdges *theEdges = femEdgesCreate(theProblem->mesh);
        for (i = 0; i < theEdges->nEdge; i++) {      
            if (theEdges->edges[i].elem[1] < 0) {       
                theProblem->dirichlet[theEdges->edges[i].node[0]] = 1; 
                theProblem->dirichlet[theEdges->edges[i].node[1]] = 1; }}
        femEdgesFree(theEdges);
    //}
    //----------------------------------------------

    switch (solverType) {
        case FEM_FULL : 
                theProblem->solver = femSolverFullCreate(theProblem->size,
                                                         theProblem->sizeLoc); break;
        case FEM_BAND : 
                band = femMeshComputeBand(theProblem->mesh);
                theProblem->solver = femSolverBandCreate(theProblem->size,
                                                         theProblem->sizeLoc,band); break;
        case FEM_ITER : 
               theProblem->solver = femSolverIterativeCreate(theProblem->size,
                                                             theProblem->sizeLoc); break;
        default : Error("Unexpected solver option"); }
    
    theProblem->soluce = malloc(sizeof(double)*theProblem->size);
    for (i = 0; i < theProblem->size; i++) 
        theProblem->soluce[i] = 0;
    
    return theProblem;
}

void femDiffusionFree_keepMesh(femDiffusionProblem *theProblem)
{
    femIntegrationFree(theProblem->rule);
    femDiscreteFree(theProblem->space);
    //femMeshFree(theProblem->mesh);
    femSolverFree(theProblem->solver);
    free(theProblem->dirichlet);
    free(theProblem->soluce);
    free(theProblem);
}

void femDiffusionCompute_motor(femDiffusionProblem *theProblem, motor *theMotor) //motor en arg
{
    femMesh *theMesh = theProblem->mesh;
    femIntegration *theRule = theProblem->rule;
    femDiscrete *theSpace = theProblem->space;
    femSolver *theSolver = theProblem->solver;
    int *number = theProblem->mesh->number;
    //double source = theProblem->sourceValue;
    double dirichlet = theProblem->dirichletValue;

    if (theSpace->n > 4) Error("Unexpected discrete space size !");    
    double Xloc[4],Yloc[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    double Uloc[4];
    int iEdge,iElem,iInteg,i,j,map[4],ctr[4];
    double **A = theSolver->local->A;
    double *Aloc = theSolver->local->A[0];
    double *Bloc = theSolver->local->B;
    
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        //le bon mu et bon source j_s
        int iDomain = theMotor->mesh->domain[iElem];
        double mu = theMotor->mu[iDomain];
        double source = theMotor->js[iDomain]; //js terme source

        for (i = 0; i < theSpace->n; i++)      Bloc[i] = 0;
        for (i = 0; i < (theSpace->n)*(theSpace->n); i++) Aloc[i] = 0;
        femDiffusionMeshLocal(theProblem,iElem,map,ctr,Xloc,Yloc,Uloc);  
        for (iInteg=0; iInteg < theRule->n; iInteg++) {    
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];  
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
            double dxdxsi = 0;
            double dxdeta = 0;
            double dydxsi = 0; 
            double dydeta = 0;
            for (i = 0; i < theSpace->n; i++) {    
                dxdxsi += Xloc[i]*dphidxsi[i];       
                dxdeta += Xloc[i]*dphideta[i];   
                dydxsi += Yloc[i]*dphidxsi[i];   
                dydeta += Yloc[i]*dphideta[i]; }
            double jac = dxdxsi * dydeta - dxdeta * dydxsi;
            for (i = 0; i < theSpace->n; i++) {    
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; }            
            for (i = 0; i < theSpace->n; i++) { 
                for(j = 0; j < theSpace->n; j++) {
                    //A_ij * 1/mu
                    A[i][j] += (1.0/ mu ) * (
                        (dphidx[i] * dphidx[j] 
                              + dphidy[i] * dphidy[j]) * jac * weight)
                    ; 
                }
            }                                                                                            
            for (i = 0; i < theSpace->n; i++) {
                Bloc[i] += phi[i] * jac * source * weight; }}
        for (i = 0; i < theSpace->n; i++) 
            if (ctr[i] == 1) femFullSystemConstrain(theSolver->local,i,dirichlet);
        femSolverAssemble(theSolver,Aloc,Bloc,Uloc,map,theSpace->n); } 
 
    double *soluce = femSolverEliminate(theSolver);
    for (i = 0; i < theProblem->size; i++)
        theProblem->soluce[i] += soluce[number[i]];
}

//Q2
double theta_node(int node_index){
    double X1 = theMesh->X[node_index];
    double Y1 = theMesh->Y[node_index];
    printf("X1:%f Y1:%f  / ",X1,Y1);
    return 2 * atan(Y1/ ( X1+sqrtf(pow(X1,2) + pow(Y1,2))));
}

//Q2
double theta_xy(double x,double y){
    return 2 * atan(y/ ( x+sqrtf(pow(x,2) + pow(y,2))));
}

int cmp_nodes(const void * a, const void * b) {
    int n1 = *(int*) a;
    int n2 = *(int*) b;

    double X1 = theMesh->X[n1]; double X2 = theMesh->X[n1];
    double Y1 = theMesh->Y[n2]; double Y2 = theMesh->Y[n2];
    double theta1 = 2 * atan(Y1/ ( X1+sqrtf(pow(X1,2) + pow(Y1,2))));
    double theta2 = 2 * atan(Y2/ ( X2+sqrtf(pow(X2,2) + pow(Y2,2))));
    //printf("n1:%d X:%f Y:%f thetha1:%f //n2:%d X:%f Y:%f theta2:%f\n",n1->index,X1,Y1,theta1,n2->index,X2,Y2,theta2);
    return theta1 - theta2;
}

//le jac/ perimetre^2 
double compute_jacSurPerim(motorMesh* theMesh, int n0,int n1,int n2){
        double* X = theMesh->X;
        double* Y = theMesh->Y;
        double localJac = (X[n1]-X[n0])*(Y[n2]-Y[n0]) - (Y[n1]-Y[n0])*(X[n2]-X[n0]);
        double edge1 = sqrt(pow(X[n1]-X[n0],2) + pow(Y[n1]-Y[n0],2));
        double edge2 = sqrt(pow(X[n1]-X[n2],2) + pow(Y[n1]-Y[n2],2));;
        double edge3 = sqrt(pow(X[n2]-X[n0],2) + pow(Y[n2]-Y[n0],2));;
        double perimeter = pow(edge1 + edge2 + edge3,3);
        return fabs(localJac) / perimeter;
}

//
// ========= Projet � r�aliser ===================
//
/*
la partie a opti c'est le solveur de ComputeMagneticPotential 
pcq c'est bcq plus long que toutes les autres fonctions ordre de plus de 10x ou 100x


*/

void motorAdaptMeshOpti(motor *theMotor, double delta)
{
    motorMesh *theMesh = theMotor->mesh;

    double x,y;
    for(int iNode = 0; iNode < theMesh->nNode; ++iNode){
        if(theMotor->movingNodes[iNode] == 1){
            //tourne partie déja présente
            x = theMesh->X[iNode]*cos(delta) - theMesh->Y[iNode]*sin(delta);
            y = theMesh->X[iNode]*sin(delta) + theMesh->Y[iNode]*cos(delta);
            theMesh->X[iNode] = x;
            theMesh->Y[iNode] = y; 
        }
    }
    theMotor->theta += delta;

    //nouveau
    
    double* X = theMesh->X;
    double* Y = theMesh->Y;
    int nElem_Airgap = theMesh->nElemDomain[11];

    int nbre_noeud_airgap = 3 * theMesh->nElemDomain[11]; //magie jsp pq autant d'elem que de noeuds
    int* statorNodes = malloc(sizeof(int)*nbre_noeud_airgap);
    uint len_statorNodes = 0;
    int* rotorNodes = malloc(sizeof(int)*nbre_noeud_airgap);
    uint len_rotorNodes = 0;

    //rempli liste statornode et rotornode avec index des noeuds correspondants
    for (size_t iElem = theMesh->nElem - nElem_Airgap; iElem < theMesh->nElem; iElem++){
        for (size_t iNode = 0; iNode < 3; iNode++){
            int node = theMesh->elem[iElem*3 + iNode];/////
            
            if(theMotor->movingNodes[node] == 1){ //air_gap partie rotor
                rotorNodes[len_rotorNodes++] = node;
            }else{          //air_gap partie stator
                statorNodes[len_statorNodes++] = node;
            }
        }
    }
    
    //tri en fonction de leur angle
    qsort(statorNodes,len_statorNodes,sizeof(int),cmp_nodes);
    qsort(rotorNodes,len_rotorNodes,sizeof(int),cmp_nodes);

    //print des triés
    printf("\n statornodes: len:%u\n",len_statorNodes);
    for (size_t i = 0; i < len_statorNodes; i++){
        printf("%d, ",statorNodes[i]);
    }printf("\n rotornodes: len:%u\n",len_rotorNodes);
    for (size_t i = 0; i < len_rotorNodes; i++){
        printf("%d, ",rotorNodes[i]);
    }fflush(stdout);
    
    printf("REMAILLAGE\n");
//REMAILLAGE NOUVEAU ELEM
    uint i_newElem = theMesh->nElem - nElem_Airgap;

    //remaillage triangle avec edge stator 
    uint i_rotor = 0;
    //pour chaque edge de statornodes
    for (size_t i_stator = 0; i_stator < len_statorNodes; i_stator++)
    {
        int i_statornext = (i_stator + 1) % len_statorNodes;
        //cherche le meilleur Jac
        double nexJac = compute_jacSurPerim(theMesh,statorNodes[i_stator],statorNodes[i_statornext],rotorNodes[i_rotor++]);
        double oldJac = 0.0;
        printf("i0:%ld i1:%d i2:%d n0:%d n1:%d n2:%d oldJac:%f nexJac=%f\n",i_stator, i_rotor, i_statornext, statorNodes[i_stator],rotorNodes[i_rotor-1],statorNodes[i_statornext],oldJac,nexJac);
        while(nexJac >= oldJac){
             oldJac = nexJac;
            nexJac = compute_jacSurPerim(theMesh,statorNodes[i_stator],statorNodes[i_statornext],rotorNodes[i_rotor++]);
            printf("i0:%ld i1:%d i2:%d n0:%d n1:%d n2:%d oldJac:%f nexJac=%f\n",i_stator, i_rotor, i_statornext, statorNodes[i_stator],rotorNodes[i_rotor-1],statorNodes[i_statornext],oldJac,nexJac);
          
        }

        //note le nouveau triangle dans elem
        theMesh->elem[3 * i_newElem] = statorNodes[i_stator];
        theMesh->elem[3 * i_newElem + 1] = rotorNodes[i_rotor - 1];
        theMesh->elem[3 * i_newElem + 2] = statorNodes[i_statornext];

        i_newElem++;
        oldJac = 0.0;
    }

    //remaillage triangle avec edge rotor 
    uint i_stator = 0;
    //pour chaque edge de statornodes
    for (size_t i_rotor = 0; i_rotor < len_rotorNodes; i_rotor++)
    {
        int i_rotornext = (i_rotor + 1) % len_rotorNodes;
        //cherche le meilleur Jac
        double nexJac = compute_jacSurPerim(theMesh,statorNodes[i_rotor],rotorNodes[i_rotornext],statorNodes[i_stator++]);
        double oldJac = 0.0;

        while(nexJac >= oldJac){
            oldJac = nexJac;
            nexJac = compute_jacSurPerim(theMesh,statorNodes[i_rotor],statorNodes[i_rotornext],rotorNodes[i_stator++]);
        }

        //note le nouveau triangle dans elem
        theMesh->elem[3 * i_newElem] = rotorNodes[i_rotor];
        theMesh->elem[3 * i_newElem + 1] = statorNodes[i_stator - 1];
        theMesh->elem[3 * i_newElem + 2] = rotorNodes[i_rotornext];

        i_newElem++;
    }
    
    

}

void motorAdaptMesh_simple(motor *theMotor, double delta)
{
    if(theMesh == NULL){
        theMesh = theMotor->mesh;
        theMotor_ = theMotor;
    }

    double x,y;
    for(int iNode = 0; iNode < theMesh->nNode; ++iNode){
        if(theMotor->movingNodes[iNode] == 1){
            //tourne partie déja présente
            x = theMesh->X[iNode]*cos(delta) - theMesh->Y[iNode]*sin(delta);
            y = theMesh->X[iNode]*sin(delta) + theMesh->Y[iNode]*cos(delta);
            theMesh->X[iNode] = x;
            theMesh->Y[iNode] = y; 
        }
    }
    theMotor->theta += delta;

    //nouveau
    
    double* X = theMesh->X;
    double* Y = theMesh->Y;
    int nElem_Airgap = theMesh->nElemDomain[11];

    int nbre_noeud_airgap = 3 * theMesh->nElemDomain[11]; //magie jsp pq autant d'elem que de noeuds
    int* statorNodes = malloc(sizeof(int)*nbre_noeud_airgap);
    uint len_statorNodes = 0;
    int* rotorNodes = malloc(sizeof(int)*nbre_noeud_airgap);
    uint len_rotorNodes = 0;

    //rempli liste statornode et rotornode avec index des noeuds correspondants
    for (size_t iElem = theMesh->nElem - nElem_Airgap; iElem < theMesh->nElem; iElem++){
        for (size_t iNode = 0; iNode < 3; iNode++){
            int node = theMesh->elem[iElem*3 + iNode];/////
            
            if(theMotor->movingNodes[node] == 1){ //air_gap partie rotor
                rotorNodes[len_rotorNodes++] = node;
            }else{          //air_gap partie stator
                statorNodes[len_statorNodes++] = node;
            }
        }
    }

    //print des triés
    /*
    printf("\n statornodes: len:%u\n",len_statorNodes);
    for (size_t i = 0; i < len_statorNodes; i++){
        printf("%d, ",statorNodes[i]);
    }printf("\n rotornodes: len:%u\n",len_rotorNodes);
    for (size_t i = 0; i < len_rotorNodes; i++){
        printf("%d, ",rotorNodes[i]);
    }fflush(stdout);*/
    
    //printf("\nREMAILLAGE\n");
//REMAILLAGE NOUVEAU ELEM
    

    //remaillage triangle avec edge stator 
    
    //pour chaque edge de statornodes
    for (size_t iElem = theMesh->nElem - nElem_Airgap; iElem < theMesh->nElem; iElem++){
        int n0 = theMesh->elem[3 * iElem];
        int n1 = theMesh->elem[3 * iElem + 1];
        int n2 = theMesh->elem[3 * iElem + 2];

        double best_jac = 0.0;
        int best_node_index = 0;

        if(theMotor_->movingNodes[n1] == 1){//edge dans le stator
            for (size_t i_rotor = 0; i_rotor < len_rotorNodes; i_rotor++)
            {
                double new_Jack = compute_jacSurPerim(theMesh,n0,rotorNodes[i_rotor],n2);
                if (new_Jack > best_jac)
                {
                    best_node_index = i_rotor;
                    best_jac = new_Jack;
                }
            }
            //note le nouveau triangle dans elem
            theMesh->elem[3 * iElem + 1] = rotorNodes[best_node_index];
        }else{
            for (size_t i_stator = 0; i_stator < len_statorNodes; i_stator++)
            {
                double new_Jack = compute_jacSurPerim(theMesh,n0,statorNodes[i_stator],n2);
                if (new_Jack > best_jac)
                {
                    best_node_index = i_stator;
                    best_jac = new_Jack;
                }
            }
            //note le nouveau triangle dans elem
            theMesh->elem[3 * iElem + 1] = statorNodes[best_node_index];
        }

        
        
        /*
        printf("i_newElem:%ld n0:%d(%f) n1:%d(%f) n2:%d(%f) \n",iElem,
        theMesh->elem[3 * iElem],theta(theMesh->elem[3 * iElem]),
        theMesh->elem[3 * iElem+1],theta(theMesh->elem[3 * iElem+1]),
        theMesh->elem[3 * iElem+2], theta(theMesh->elem[3 * iElem+2])
        );*/
    }
}

void motorAdaptMesh(motor *theMotor, double delta){
    motorAdaptMesh_simple(theMotor, delta);
}

/*
400   : Couple : -0.011850554805193 d : 0.000333333203504
838   : Couple : -0.027798564474380 d : 0.000333333723617
1667  : Couple : -0.038487112747265 d : 0.000333333037516
4424  : Couple : -0.044461436111684 d : 0.000333332935723
14608 : Couple : -0.045267709217747 d : 0.000333332667550
*/

double motorComputeCouple(motor *theMotor)
{
    double jac;
    double I = 0;
    double xLoc[3], x[3];
    double yLoc[3], y[3];
    double r[3];
    int iElem, i, j;

    static const double xsi[3]     = { 0.166666666666667, 0.666666666666667, 0.166666666666667};
    static const double eta[3]     = { 0.166666666666667, 0.166666666666667, 0.666666666666667};
    static const double weight[3]  = { 0.166666666666667, 0.166666666666667, 0.166666666666667};
    int n = 3;
    
    double dphidxsi[3] = {-1.0, 1.0, 0.0};
    double dphideta[3] = {-1.0, 0.0, 1.0};
    
    int first_rotorgap = theMesh->nElem - theMesh->nElemDomain[10] - theMesh->nElemDomain[11];
    int last_rotorgap = theMesh->nElem - theMesh->nElemDomain[11];

    double max_r = -1;
    double min_r = 1000;

    
    for(iElem = first_rotorgap; iElem < last_rotorgap; iElem++) {
        //x, y, et rayon 
        for(i = 0; i < 3; i++) {
            x[i] = theMesh->X[theMesh->elem[3*iElem+i]] ;
            y[i] = theMesh->Y[theMesh->elem[3*iElem+i]] ; 
            double rayon = sqrt(pow(x[i],2) + pow(y[i],2));
            if(rayon>max_r){ max_r = rayon;
            }else if(rayon<min_r){ min_r = rayon; }
            r[i] = rayon;
        }
        //xLoc yLoc
        for(i = 0; i < n; i++) {
            xLoc[i] = x[0] * (1- xsi[i] - eta[i]) + x[1] * xsi[i] + x[2] * eta[i];
            yLoc[i] = y[0] * (1- xsi[i] - eta[i]) + y[1] * xsi[i] + y[2] * eta[i];
        }
        //calcule deriveé de x,y en fct de xsi,eta
        double dxdxsi = 0;
        double dxdeta = 0;
        double dydxsi = 0;
        double dydeta = 0;
        for (i = 0; i < 3; i++) {  
            dxdxsi += x[i]*dphidxsi[i];       
            dxdeta += x[i]*dphideta[i];   
            dydxsi += y[i]*dphidxsi[i];   
            dydeta += y[i]*dphideta[i]; 
        }
        
        jac = dxdxsi * dydeta - dxdeta * dydxsi;

        //calcul de da/dx et da/dy
        double dphidx = 0.0;
        double dphidy = 0.0;
        double dadx = 0.0;
        double dady = 0.0;
        for (i = 0; i < 3; i++) {    
            dphidx = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
            dphidy = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; 
            double a = theMotor->a[theMesh->elem[3*iElem+i]];
            dadx += dphidx*a;
            dady += dphidy*a;
        }

        //calcul de l'intégrale
        for (size_t i_node = 0; i_node < 3; i_node++){
            double i_added = weight[i_node] * (pow(10,7)/r[i_node]) * (dadx * x[i_node] + dady * y[i_node]) * (dady * x[i_node] - dadx * y[i_node]);
            I = i_added;
            //printf("/i_added:%f (x,y):(%f,%f) r:%f \n",i_added,x[i_node],y[i_node],r[i_node]);
        }

        //printf("        ielem:%d dadx:%f dady:%f I:%f\n\n",iElem,dadx,dady,I);

        
    }
    double d = max_r - min_r;//a changer
    double mu_0 = 4*M_PI;
    double cst =  theMotor_-> L / (d * mu_0);
    
    I = - I * cst;
    printf("d=%f cst:%f I=%f \n",d,cst,I);
    return 0.0;
}

double motorComputeCouple_dru(motor *theMotor)
{
    //return 0.0;

    double jac;
    double I = 0;
    double xLoc[3], x[3];
    double yLoc[3], y[3];
    double r[3];
    int iElem, i, j;

    static const double xsi[3]     = { 0.166666666666667, 0.666666666666667, 0.166666666666667};
    static const double eta[3]     = { 0.166666666666667, 0.166666666666667, 0.666666666666667};
    static const double weight[3]  = { 0.166666666666667, 0.166666666666667, 0.166666666666667};
    int n = 3;
    /*
    double x2xsi[3] =  {0.0, 1.0, 0.0};
    double x2eta[3] =  {0.0, 0.0, 1.0};
    */
    
    double dphidxsi[3] = {-1.0, 1.0, 0.0};
    double dphideta[3] = {-1.0, 0.0, 1.0};
    
    int first_rotorgap = theMesh->nElem - theMesh->nElemDomain[10] - theMesh->nElemDomain[11];
    int last_rotorgap = theMesh->nElem - theMesh->nElemDomain[11];

    double max_r = -1;
    double min_r = 1000;

    double I_tim = 0.0;
    
    for(iElem = first_rotorgap; iElem < last_rotorgap; iElem++) {
        //x et y position 
        for(i = 0; i < 3; i++) {
            x[i] = theMesh->X[theMesh->elem[3*iElem+i]] ;
            y[i] = theMesh->Y[theMesh->elem[3*iElem+i]] ; 
            double rayon = sqrt(pow(x[i],2) + pow(y[i],2));
            if(rayon>max_r){ max_r = rayon;
            }else if(rayon<min_r){ min_r = rayon; }
            r[i] = rayon;
            printf("rayon:%f    ",rayon);
        }
        //xLoc yLoc
        for(i = 0; i < n; i++) {
            xLoc[i] = x[0] * (1- xsi[i] - eta[i]) + x[1] * xsi[i] + x[2] * eta[i];
            yLoc[i] = y[0] * (1- xsi[i] - eta[i]) + y[1] * xsi[i] + y[2] * eta[i];
        }
        //calcule deriveé de x,y en fct de xsi,eta
        double dxdxsi = 0;
        double dxdeta = 0;
        double dydxsi = 0;
        double dydeta = 0;
        for (i = 0; i < 3; i++) {  
            dxdxsi += x[i]*dphidxsi[i];       
            dxdeta += x[i]*dphideta[i];   
            dydxsi += y[i]*dphidxsi[i];   
            dydeta += y[i]*dphideta[i]; 
        }
        
        jac = dxdxsi * dydeta - dxdeta * dydxsi;

        //calcul de da/dx et da/dy
        double dphidx = 0.0;
        double dphidy = 0.0;
        double dadx = 0.0;
        double dady = 0.0;
        for (i = 0; i < 3; i++) {    
            dphidx = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
            dphidy = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; 
            double a = theMotor->a[theMesh->elem[3*iElem+i]];
            dadx += dphidx*a;
            dady += dphidy*a;
            //printf("index:%d a:%f\n",theMesh->elem[3*iElem+i],a);
            
        }

        //calcul de l'intégrale
        double add_i_amine = 0.0;
        for (size_t i_node = 0; i_node < 3; i_node++){
            add_i_amine += weight[i_node] * (pow(10,7)/r[i_node]) * (dadx * x[i_node] + dady * y[i_node]) * (dady * x[i_node] - dadx * y[i_node]);
            printf("add_i_amine:%f (x,y):(%f,%f) r:%f / ",add_i_amine,x[i_node],y[i_node],r[i_node]);
        }
        printf("\n");

        //calcul du terme a ajouter a l'intégrale
        double x_moyen = (x[0] + x[1] + x[2]) / 3;
        double y_moyen = (y[0] + y[1] + y[2]) / 3;
        
        double r_moyen = sqrt(pow(x_moyen,2) + pow(y_moyen,2));

        double angle = theta_xy(x_moyen,y_moyen);

        double add_i_trigo = pow(10,7)* (dadx * cos(angle) + dady * sin(angle)) * (dady * r_moyen * cos(angle) - dadx * r_moyen * sin(angle));
        double add_i_xy = (pow(10,7)/r_moyen) * (dadx * x_moyen + dady * y_moyen) * (dady * x_moyen - dadx * y_moyen);
        
        I_tim += add_i_amine;
        printf("ielem:%d x_moyen:%f y_moyen:%f dadx:%f dady:%f angle:%f cos(a):%f sin(a):%f r_moyen:%f +amne:%f +trigo:%f +xy:%f I_tim:%f\n",
            iElem,x_moyen,y_moyen,dadx,dady,angle,cos(angle),sin(angle),
            r_moyen,add_i_amine,add_i_trigo,add_i_xy,I_tim
        );

        
    }
    double d = max_r - min_r;//a changer
    double mu_0 = 4*M_PI;
    double cst =  theMotor_-> L / (d * mu_0);
    
    I = - I * cst * pow(10,7);
    I_tim = - I_tim * cst;
    printf("d=%f cst:%f I=%f I_tim=%f\n",d,cst,I,I_tim);
    return 0.0;
}



//20 <= theta < 50 et 50 <= theta < 80
//quel que soit langle de theta on le ramene au premier quadrant
//activer les bobines tjrs par paires A+ avec A- etc
void motorComputeCurrent(motor *theMotor)
{
    double theta = theMotor->theta;
    double prem_quadr = M_PI / 2;
    //printf("theta:%f prem_quadr:%f M_PI:%f",theta,prem_quadr,M_PI);
    while(theta > 0){
        theta -= prem_quadr;
    }
    theta += prem_quadr;
    printf(" prem _quadrant:%f",theta);

    if(20.0/180 * M_PI <= theta && theta < 50.0/180 * M_PI){ //coils A (1,2)
        printf(" -> 20<=theta<50");
        theMotor->js[1] = 8.8464 * pow(10,5);
        theMotor->js[2] = - 8.8464 * pow(10,5);
        theMotor->js[3] = 0.0;
        theMotor->js[4] = 0.0;
        theMotor->js[5] = 0.0;
        theMotor->js[6] = 0.0;
    }else if(50.0/180 * M_PI <= theta && theta < 80.0/180 * M_PI){ //coils B (3,4)
        printf(" -> 50<=theta<80");
        theMotor->js[1] = 0.0;
        theMotor->js[2] = 0.0;
        theMotor->js[3] = 8.8464 * pow(10,5);
        theMotor->js[4] = - 8.8464 * pow(10,5);
        theMotor->js[5] = 0.0;
        theMotor->js[6] = 0.0;

    }else{                                                   //coils C (5,6)
        printf(" ->  80<=theta<20");
        theMotor->js[1] = 0.0;
        theMotor->js[2] = 0.0;
        theMotor->js[3] = 0.0;
        theMotor->js[4] = 0.0;
        theMotor->js[3] = 8.8464 * pow(10,5);
        theMotor->js[4] = - 8.8464 * pow(10,5);
    }
    printf("\n");
    return;
}

void motorComputeMagneticPotential(motor *theMotor)
{
    if(theMesh == NULL){
        theMesh = theMotor->mesh;
        theMotor_ = theMotor;
    }


    femMesh *thefemMesh = motorMesh_to_femMesh(theMotor->mesh);
    
    femDiffusionProblem* theProblem = femDiffusionCreate_mesh(thefemMesh,FEM_BAND,FEM_XNUM);
    

    femDiffusionCompute_motor(theProblem,theMotor);

    for (size_t i = 0; i < theMesh->nNode; i++)
    {
        theMotor->a[i] = theProblem->soluce[i];
    }
    
    

    femDiffusionFree_keepMesh(theProblem);

    
    return;
}

//
// ========= Projet � r�aliser ===================
//







