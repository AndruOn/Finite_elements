#include "fem.h"
#include "motor.h"

//variables pour performance
double offset = 0.2;

//Variable utile
motorMesh* theMesh;
motor* theMotor_;
femDiffusionProblem* theProblem;

//variable utile pour la récolte de statistiques
char* stats_filename = "../data/mesh4424.csv";
FILE* file;

//time functiuns
#include <time.h>
clock_t start, end;
double cpu_time_used;

//variable de print
int DEBUG_PRINT = 0;
int WRITE_STATS = 0;
    
/*
Ecrit dans le fichier de sortie stats_filename pour l'écriture de données
*/
int create_stats(){
    if(file == NULL){
        file = fopen(stats_filename,"w+");
        if(file == NULL){
            printf("ERROR: coulnd fopen(stats_filename%s","");
            return -1;
        }
        //fprintf(file,"cputime,time,couple,theta,omega,coil,offset\n");
        theMotor_->omega = 0.0;
        if(DEBUG_PRINT) printf("FILE STATS CREATED\n");
    }
    return 0;
}

//
// ========= Fonctions intermédiaires ===================
//

//Q1--------------------------------------------------------------------------

/*
Fait un cast "propre" de motorMesh en femMesh
Retourne une structure femMesh représentant le même mesh que dans motorMesh donné en argument
*/
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

/*
femDiffusionFree modifiée pour ne pas free le mesh à chaque étape
Cette fonction crée une structure femDiffusionProblem pour le calcul du potentiel magnétique 
*/
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
    for (i = 0; i < theProblem->size; i++) {
        theProblem->dirichlet[i] = 0;
    }
    femEdges *theEdges = femEdgesCreate(theProblem->mesh);
    for (i = 0; i < theEdges->nEdge; i++) {      
        if (theEdges->edges[i].elem[1] < 0) {       
            theProblem->dirichlet[theEdges->edges[i].node[0]] = 1; 
            theProblem->dirichlet[theEdges->edges[i].node[1]] = 1; 
        }
    }
    femEdgesFree(theEdges);

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

/*
femDiffusionFree modifiée pour ne pas free le mesh à chaque étape
Cette fonction free la struture femDifussionProblem
*/
void femDiffusionFree_keepMesh(femDiffusionProblem *theProblem)
{
    femIntegrationFree(theProblem->rule);
    femDiscreteFree(theProblem->space);
    femSolverFree(theProblem->solver);
    free(theProblem->dirichlet);
    free(theProblem->soluce);
    free(theProblem->mesh->number);
    free(theProblem->mesh);
    free(theProblem);
}

/*
Nettoie la struture theProblem pour la prochaine itération
free et remalloc le solveur ainsi que remettre solucs à zéro
*/
void clean_theProblem(femDiffusionProblem* theProblem,femSolverType solverType){
    for (int i = 0; i < theProblem->size; i++) 
        theProblem->soluce[i] = 0;

    femSolverFree(theProblem->solver);
    int band;
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
}

/*
femDiffusionCompute modifiée pour calculer le potentiel magnétique 
de la structure motor donné en argument
*/
void femDiffusionCompute_motor(femDiffusionProblem *theProblem, motor *theMotor) 
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
        theMotor->a[i] = soluce[number[i]];
}


//Q2--------------------------------------------------------------------------

/*
Retourne l'angle en radian (strictement positif:[0:2*M_PI])
d'un point dont on donne la position (x,y) en argument
*/
double theta_xy(double x,double y){
    return 2 * atan(y/ ( x+sqrtf(pow(x,2) + pow(y,2))));
}

//Retourne le jac/perimetre^3 d'un triangle en prenant comme argument les 3 noueds du triangle et le mesh qui contient ce triangle
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


//Q4--------------------------------------------------------------------------

/*
Allume le courrant dans la paire de bobines données en argument, l'éteint dans les autres paires de bobines
*/
void turnOn_current_in_coil(char c, motor* theMotor){
    switch (c)
    {
    case 'A':
        theMotor->js[1] = 8.8464 * pow(10,5);
        theMotor->js[2] = - 8.8464 * pow(10,5);
        theMotor->js[3] = 0.0;
        theMotor->js[4] = 0.0;
        theMotor->js[5] = 0.0;
        theMotor->js[6] = 0.0;
        break;
    case 'B':
        theMotor->js[1] = 0.0;
        theMotor->js[2] = 0.0;
        theMotor->js[3] = 8.8464 * pow(10,5);
        theMotor->js[4] = - 8.8464 * pow(10,5);
        theMotor->js[5] = 0.0;
        theMotor->js[6] = 0.0;

        break;
    case 'C':
        theMotor->js[1] = 0.0;
        theMotor->js[2] = 0.0;
        theMotor->js[3] = 0.0;
        theMotor->js[4] = 0.0;
        theMotor->js[5] = 8.8464 * pow(10,5);
        theMotor->js[6] = - 8.8464 * pow(10,5);
        break;
    
    default:
        printf("ERROR: Only Coils A, B and C exist \n");
        break;
    }
}

//
// ========= Fonctions appelées dans main lors d'une itération ===================
//

/*
Q2
Modifie le mesh de theMotor en faisant tourner les noueds du
rotor de delta radians et remaille les triangle du domaine Air_gap 
*/
void motorAdaptMesh(motor *theMotor, double delta)
{
    if(theMesh == NULL){
        theMesh = theMotor->mesh;
        theMotor_ = theMotor;
    }

    //Fais tourner tout les noueds du rotor
    double x,y;
    for(int iNode = 0; iNode < theMesh->nNode; ++iNode){
        if(theMotor->movingNodes[iNode] == 1){
            x = theMesh->X[iNode]*cos(delta) - theMesh->Y[iNode]*sin(delta);
            y = theMesh->X[iNode]*sin(delta) + theMesh->Y[iNode]*cos(delta);
            theMesh->X[iNode] = x;
            theMesh->Y[iNode] = y; 
        }
    }
    theMotor->theta += delta;

    
    double* X = theMesh->X;
    double* Y = theMesh->Y;
    int nElem_Airgap = theMesh->nElemDomain[11];

    int nbre_noeud_airgap = 3 * theMesh->nElemDomain[11]; //autant d'elem que de noeuds et 3* car 3 noueds par trinagle
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
    
    //pour chaque elem change le noeud du milieu pour tester le triangle avec le meilleur jac/perim^3
    for (size_t iElem = theMesh->nElem - nElem_Airgap; iElem < theMesh->nElem; iElem++){
        int n0 = theMesh->elem[3 * iElem];
        int n1 = theMesh->elem[3 * iElem + 1];
        int n2 = theMesh->elem[3 * iElem + 2];

        double best_jac = 0.0;
        int best_node_index = 0;

        if(theMotor_->movingNodes[n1] == 1){//edge(lesdeux noueds qui restent connectés) dans le stator
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
        }else{//edge(lesdeux noueds qui restent connectés) dans le rotor
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
    }
    free(statorNodes);
    free(rotorNodes);
}



/*
Q3
Retourne le couple appliquée au rotor par le potentiel calculé sur chaque noeud du motorMesh
le couple est calculé avec une règle d'intégration de Hammer 3 points
*/
double motorComputeCouple(motor *theMotor){
    double jac;
    double I = 0.0;
    double x[3], y[3];
    int nodes[3];
    double xLoc, yLoc;
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
            nodes[i] = theMesh->elem[3*iElem+i];
            x[i] = theMesh->X[nodes[i]] ;
            y[i] = theMesh->Y[nodes[i]] ; 
            double rayon = sqrt(pow(x[i],2) + pow(y[i],2));
            if(rayon>max_r){ max_r = rayon;
            }else if(rayon<min_r){ min_r = rayon; }
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
        
        jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);

        //calcul de da/dx et da/dy
        double dphidx = 0.0;
        double dphidy = 0.0;
        double dadx = 0.0;
        double dady = 0.0;
        for (i = 0; i < 3; i++) {    
            dphidx = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
            dphidy = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; 
            double a = theMotor->a[nodes[i]];
            dadx += dphidx*a;
            dady += dphidy*a;
        }

        //calcul de l'intégrale
        for (size_t i_node = 0; i_node < 3; i_node++){
            xLoc = x[0]*(1.0 - xsi[i_node] - eta[i_node]) + x[1]*xsi[i_node] + x[2]*eta[i_node];
            yLoc= y[0]*(1.0 - xsi[i_node] - eta[i_node]) + y[1]*xsi[i_node] + y[2]*eta[i_node];
            double r = sqrt(pow(xLoc,2) + pow(yLoc,2));;
            double dadtheta_dadr = (dadx * xLoc + dady * yLoc) * (dady * xLoc /r - dadx * yLoc/r);
            double i_added = jac * weight[i_node] * dadtheta_dadr;
            I += i_added;
        }
    }

    double d = max_r - min_r;//a changer
    double mu_0 = 4*M_PI*pow(10,-7);
    double cst =  theMotor-> L / (d * 4*M_PI*1e-7);
    
    I = - I * cst;

    if(DEBUG_PRINT) printf("3: Couple=%f d=%f rmax:%f rmin:%f\n",I,d,max_r,min_r);
    if(WRITE_STATS) fprintf(file,"%f,%0.10f,%0.6f,%0.6f,",theMotor_->time,I,theMotor_->theta,theMotor_->omega);
    return I;
}



/*
Q4
Allume la paire de bobines adéquate en fonction de l'angle(theta) de la structure theMotor
J'ajoute un offset à l'angle et le ramène au premier quadrant avant de vérifier
dans quel tiers du premier quadrant ce angle ce retoruve pour allumer 
la paire de bobines correspondante
*/
void motorComputeCurrent(motor *theMotor)
{   
    double theta = theMotor->theta;
    if(DEBUG_PRINT) printf("4:theta:%f ",theta);
    //0.172 gauthier trop loin
    //20.0/180 * M_PI andru
    

    theta+=offset;

    double prem_quadr = M_PI / 2;
    if(theta>0){
        while(theta > 0){
            theta -= prem_quadr;
        }
        theta += prem_quadr;
    }else{
        while(theta <= 0){
            theta += prem_quadr;
        }
    }
    if(DEBUG_PRINT) printf("offset:%f theta mod pi/2:%f ", offset, theta);
   
    if(0.0 <= theta && theta < M_PI/6){ //coils B 
        char coil = 'B';
        if(DEBUG_PRINT) printf(" -> %f<=theta<%f Coil %c\n",0.0,M_PI/6,coil);
        if(WRITE_STATS) fprintf(file,"%c,%0.5f\n",coil,offset);
        turnOn_current_in_coil(coil, theMotor);

    }else if(M_PI/6 <= theta && theta < 2 * M_PI/6){//coils A
        char coil = 'A';
        if(DEBUG_PRINT) printf(" -> %f<=theta<%f Coil %c\n", M_PI/6, 2*M_PI/6,coil);
        if(WRITE_STATS) fprintf(file,"%c,%0.5f\n",coil,offset);
        turnOn_current_in_coil(coil, theMotor);

    }else if(2 * M_PI/6 <= theta && theta < 3 * M_PI/6){ //coils C 
        char coil = 'C';
        if(DEBUG_PRINT) printf(" -> %f<=theta<%f Coil %c\n", 2*M_PI/6, 3*M_PI/6,coil);
        if(WRITE_STATS) fprintf(file,"%c,%0.5f\n",coil,offset);
        turnOn_current_in_coil(coil, theMotor);
    }else{                                                   
        if(DEBUG_PRINT) printf("ERROR: pas dans le premier quadrant---------------------\n");
        if(WRITE_STATS) fprintf(file,"ERROR,%0.5f\n",offset);
    }
    
    return;
}




/*
Q1
Calcule le potentiel magnétique sur le mesh de theMotor
Remplit Motor->a avec le potentiel magnétique
*/
void motorComputeMagneticPotential(motor *theMotor)
{
    if(DEBUG_PRINT || WRITE_STATS) start = clock();

    if(theMesh == NULL){
        theMesh = theMotor->mesh;
        theMotor_ = theMotor;
    }
    if(WRITE_STATS) create_stats();

    /*
    if(theProblem == NULL){
        femMesh *thefemMesh = motorMesh_to_femMesh(theMotor->mesh);
        theProblem = femDiffusionCreate_mesh(thefemMesh,FEM_BAND,FEM_XNUM);
    }else{
        clean_theProblem(theProblem,FEM_BAND);
    }*/


    femMesh *thefemMesh = motorMesh_to_femMesh(theMotor->mesh);
    theProblem = femDiffusionCreate_mesh(thefemMesh,FEM_BAND,FEM_XNUM);
    
    
    femDiffusionCompute_motor(theProblem,theMotor);
    femDiffusionFree_keepMesh(theProblem);
    

    if(DEBUG_PRINT){
        end = clock();
        cpu_time_used = (double) (end - start) / CLOCKS_PER_SEC;
        printf("1: time_magnetcpotential:%f\n",cpu_time_used);
    }
    if(WRITE_STATS){
        end = clock();
        cpu_time_used = (double) (end - start) / CLOCKS_PER_SEC;
        fprintf(file,"%0.7f,",cpu_time_used);
    }
    
    return;
}


/*
Free la structure Motor donnée en argument
*/
void motorFree(motor *theMotor)
{
    //femDiffusionFree_keepMesh(theProblem);
    free(theMotor->mu);
    free(theMotor->js);
    free(theMotor->movingNodes);
    free(theMotor->a);
    free(theMotor);
    return;
}
