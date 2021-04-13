# include "fem.h"


femProblem *advdiffNew(double epsilon, double beta, double zeta, int size, int sizePlot)
{
    femProblem *myProblem = malloc(sizeof(femProblem));
    myProblem->epsilon = epsilon;
    myProblem->beta    = beta;
    myProblem->zeta    = zeta;

	myProblem->size = size;    
    myProblem->X = malloc(sizeof(double) * size); 
    myProblem->U = malloc(sizeof(double) * size);
    
    myProblem->sizePlot = sizePlot;    
    myProblem->x = malloc(sizeof(double) * sizePlot); 
    myProblem->u = malloc(sizeof(double) * sizePlot); 
    myProblem->uregime = malloc(sizeof(double) * sizePlot); 
 
    advdiffReset(myProblem);


    return myProblem;
}

void advdiffReset(femProblem *myProblem)
{
    double h;
    int i;

    h = 1.0/(myProblem->size-1);
	  for (i=0; i < myProblem->size; i++) {
        myProblem->X[i] = i*h; 
        myProblem->U[i] = 0.0; }
    myProblem->U[0] = 1.0;


    h = 1.0/(myProblem->sizePlot-1);
	  for (i=0; i < myProblem->sizePlot; i++) {
        myProblem->x[i] = i*h; 
        myProblem->u[i] = 0.0;        
        myProblem->uregime[i] = 0.0; }
}

void advdiffFree(femProblem *myProblem)
{
    free(myProblem->X);
    free(myProblem->U);
    free(myProblem->x);
    free(myProblem->u);
    free(myProblem->uregime);
    free(myProblem);
}

void advdiffSteadySolution(femProblem *myProblem)
{
    int  size = myProblem->sizePlot, i;
    double *u = myProblem->uregime;
    double *x = myProblem->x;
    double beta    = myProblem->beta;
    double epsilon = myProblem->epsilon;
    double pe  = beta/epsilon; 
  
    
    if (beta != 0.0) 
       for (i=0; i < size; i++) 
          u[i] = (exp(pe) - exp(pe*x[i]) )/(exp(pe) - 1.0); 
    else 
       for (i=0; i < size; i++) 
          u[i] = 1 - x[i];          

}

# ifndef NOANALYTIC

void advdiffSolution(femProblem *myProblem, double time)
{
    int  size = myProblem->sizePlot, i,j;
    double* u = myProblem->u;
    double* x = myProblem->x;
    double beta    = myProblem->beta;
    double epsilon = myProblem->epsilon;
    double zeta = myProblem->zeta ;
    double beta_sqr = beta * beta;

    double Pe = beta/epsilon;
    double err = 1000.0;

    //Calcule sol de régime et pose u à sol régime
    advdiffSteadySolution(myProblem);
    for (i=0; i < size; i++){
        u[i] = myProblem->uregime[i];
    }

    double k,coeff,alpha;

    //U_sol = U_regime + U_transi
    for (j=1; j < 1000 && err > 1e-15 ; j++) {
        k = M_PI * j;
        coeff = ((8*k)/(Pe*Pe+4*k*k));
        alpha = beta_sqr/(4*epsilon) + k*k*epsilon;
        err = coeff *  exp(Pe/2 - alpha*time);
        for (i=0; i < size; i++) {
            u[i] = u[i] - coeff * sin(k*x[i]) * exp( -alpha*time + Pe*x[i]/2.0 ); 
        }
    }
 
  /*
  //  Solution diffusion pure :-)

    for (i=0; i < size; i++) 
        u[i] = 1 - x[i];   
    for (j=1; j < 500; j++) {
        double k = M_PI * j; 
        double alpha = k*k*epsilon;
        double coeff = 2/k;
        for (i=0; i < size; i++) {
                u[i] = u[i] - coeff * sin(k*x[i]) * exp(- alpha*time); }}
   */
    
}

# endif
# ifndef NOTIMESTEP

double advdiffComputeTimeStep(femProblem *myProblem)
{
    int size = myProblem->size;
    double beta    = myProblem->beta;
    double epsilon = myProblem->epsilon;
    double zeta = myProblem->zeta ;

    double h = 1.0/(size-1);
    double dt = fmin( (zeta*h*beta + 2*epsilon)/(beta*beta), (h*h)/(zeta*h*beta + 2*epsilon));

    //coef -> rapport entre zones de stabilite dans plan complexe
    return 1.35 * dt; 
}

# endif
# ifndef NORIGHTHANDSIDE

void advdiffComputeRightHandSide(femProblem *myProblem, double *U, double *F)
{
    int size = myProblem->size;
    double beta    = myProblem->beta;
    double epsilon = myProblem->epsilon;
    double zeta = myProblem->zeta ;
    int i;

    double h = 1.0 / (size-1);
    F[0] = 0;
    F[size-1] = 0.0;

    for (i=1; i < size-1; i++) {
        F[i] =  - ( (1-zeta) * beta * (U[i+1] - U[i-1]) / (2*h) )
            - (zeta * beta * (U[i] - U[i-1]) / h)
            + (epsilon * (U[i+1] - 2*U[i] + U[i-1]) / (h*h));
    }
    
}

# endif
# ifndef NOEULER

void advdiffUpdateEuler(femProblem *myProblem, double dt)
{
    int size = myProblem->size;

    double F[size];
    advdiffComputeRightHandSide(myProblem, myProblem->U, F);
    for (size_t i = 0; i < size; i++){
        myProblem->U[i] += dt * F[i];
    }   
}

# endif
# ifndef NORK

void advdiffUpdateRungeKutta(femProblem *myProblem, double dt)
{
    int size = myProblem->size;
    double Uold[size];

    //Keep old U
    double* U = myProblem->U;
    for (int i = 0; i < size; i++){
        Uold[i] = U[i];
    }
    
    //K1
    double K1[size];
    advdiffComputeRightHandSide(myProblem, U, K1);
    //K2
    double U2[size];
    double K2[size];
    for (size_t i = 0; i < size; i++){
        U2[i] = Uold[i] + dt/2 *K1[i];
    }
    advdiffComputeRightHandSide(myProblem, U2, K2);
    //K3
    double U3[size];
    double K3[size];
    for (size_t i = 0; i < size; i++){
        U3[i] = Uold[i] + dt/2*K2[i];
    }
    advdiffComputeRightHandSide(myProblem, U3, K3);
    //K4
    double U4[size];
    double K4[size];
    for (size_t i = 0; i < size; i++){
        U4[i] = Uold[i] + dt*K3[i];
    }
    advdiffComputeRightHandSide(myProblem, U4, K4);
    
    //U^(n+1)
    for (size_t i = 0; i < size; i++){
        U[i] += dt/6 * ( K1[i] + 2*K2[i] + 2*K3[i] + K4[i]);
    }
}

# endif

