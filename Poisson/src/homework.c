
#include"fem.h"



# ifndef NOPOISSONCREATE

femPoissonProblem *femPoissonCreate(const char *filename)
{
    femPoissonProblem *theProblem = malloc(sizeof(femPoissonProblem));
    theProblem->mesh  = femMeshRead(filename);           
    theProblem->edges = femEdgesCreate(theProblem->mesh);  
    if (theProblem->mesh->nLocalNode == 4) {
        theProblem->space = femDiscreteCreate(4,FEM_QUAD);
        theProblem->rule = femIntegrationCreate(4,FEM_QUAD); }
    else if (theProblem->mesh->nLocalNode == 3) {
        theProblem->space = femDiscreteCreate(3,FEM_TRIANGLE);
        theProblem->rule = femIntegrationCreate(3,FEM_TRIANGLE); }
    theProblem->system = femFullSystemCreate(theProblem->mesh->nNode);
    return theProblem;
}

# endif
# ifndef NOPOISSONFREE

void femPoissonFree(femPoissonProblem *theProblem)
{
    femFullSystemFree(theProblem->system);
    femIntegrationFree(theProblem->rule);
    femDiscreteFree(theProblem->space);
    femEdgesFree(theProblem->edges);
    femMeshFree(theProblem->mesh);
    free(theProblem);
}
    

# endif
# ifndef NOMESHLOCAL


void femMeshLocal(const femMesh *theMesh, const int i, int *map, double *x, double *y)
{
    int n = theMesh->nLocalNode;
    for (size_t j = 0; j < n; j++)
    {
        map[j] = theMesh->elem[i*n+j];
        x[j] = theMesh->X[map[j]];
        y[j] = theMesh->Y[map[j]];
    }
}

# endif
# ifndef NOPOISSONSOLVE


void femPoissonSolve(femPoissonProblem *theProblem)
{
  
    femMesh* mesh = theProblem->mesh;
    int* elem = mesh->elem;
    int nNode = mesh->nLocalNode;
    femIntegration* rule = theProblem->rule;
    femDiscrete* space = theProblem->space;
    double* B = theProblem->system->B;
    double** A = theProblem->system->A;
    femEdge* edges = theProblem->edges->edges;
    femFullSystem* system = theProblem->system;

    double jacobian = 0.0;
    double frontier_value = 0.0;
    double x[nNode];
    double y[nNode];
    int map[nNode];

    double phi[nNode];
    double dphi_dx[nNode];
    double dphi_dy[nNode];
    double dphi_dpsy[nNode];
    double dphi_deta[nNode];
    double dx_dpsy;
    double dx_deta;
    double dy_dpsy;
    double dy_deta;
    double xsi;
    double eta;
    double weight;

    //assemble
    for (size_t i_elem = 0; i_elem < mesh->nElem; i_elem++){
        femMeshLocal(mesh,i_elem,map,x,y);
        for (size_t i_integr = 0; i_integr < rule->n; i_integr++){
            xsi = rule->xsi[i_integr];
            eta = rule->eta[i_integr];
            weight = rule->weight[i_integr];

            femDiscretePhi2(space, xsi, eta, phi);
            femDiscreteDphi2(space, xsi, eta, dphi_dpsy, dphi_deta);
            //femDiscretePrint(space);

            dx_dpsy = 0;
            dx_deta = 0;
            dy_dpsy = 0;
            dy_deta = 0;
            //p42
            for (int i = 0;i< space->n; i++) {
                dx_dpsy += x[i] * dphi_dpsy[i];
                dx_deta += x[i] * dphi_deta[i];
                dy_dpsy += y[i] * dphi_dpsy[i];
                dy_deta += y[i] * dphi_deta[i];
            }
            jacobian = fabs(dx_dpsy * dy_deta - dx_deta * dy_dpsy);
            //p45
            for (int i = 0; i< space->n; i++) {
                dphi_dx[i] = dphi_dpsy[i] * dy_deta - dphi_deta[i] * dy_dpsy;
                dphi_dy[i] = dphi_deta[i] * dx_dpsy - dphi_dpsy[i] * dx_deta;
            }
            for (int l = 0; l < space->n; l++) {
                B[map[l]] += phi[l] * jacobian * weight;
                for (int c = 0; c < space->n; c++) {
                    A[map[l]][map[c]] += (dphi_dx[l] * dphi_dx[c] + dphi_dy[l] * dphi_dy[c]) / jacobian * weight;
                }
            }
        }
    }

    //appliquer les contraintes
    for (size_t i_elem = 0; i_elem < theProblem->edges->nEdge; i_elem++){
        femEdge edge = edges[i_elem];
        if(edge.elem[1] == -1){
            femFullSystemConstrain(system,edge.node[0],frontier_value);
            femFullSystemConstrain(system,edge.node[1],frontier_value);
        }
    }
    
    //resoudrele système
    femFullSystemEliminate(system);
}



# endif
