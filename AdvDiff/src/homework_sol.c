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
    double *u = myProblem->u;
    double *uregime = myProblem->uregime;
    double *x = myProblem->x;
    double beta    = myProblem->beta;
    double epsilon = myProblem->epsilon;
    double pe  = beta/epsilon;
    double test = 1.0;

    advdiffSteadySolution(myProblem);
    for (i=0; i < size; i++)
        u[i] = uregime[i];
    for (j=1; j < 500 && test > 1e-9 ; j++) {
        double k = M_PI * j;
        double alpha = beta*beta/(4*epsilon) + k*k*epsilon;
        double coeff = ((8*k)/(pe*pe+4*k*k));
        test = coeff *  exp(pe/2 - alpha*time);
        for (i=0; i < size; i++) {
            u[i] = u[i] - coeff * sin(k*x[i]) * exp(pe*x[i]/2.0 - alpha*time); }}

    //
    // Solution diffusion pure :-)
    //
    //   for (j=1; j < 100; j++) {
    //     double k = M_PI * j;
    //     for (i=0; i < size; i++) {
    //  		U[i] = U[i] - (2.0/k) * sin(k*X[i]) * exp(- epsilon*k*k*time); }}
    //

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

    return 1.35 * dt;   // Aussi simple que cela : eh oui :-)
    // C'est juste le rapport entre les deux zones de stabilitÃ© dans le plan complexe !
}

# endif
# ifndef NORIGHTHANDSIDE


void advdiffComputeRightHandSide(femProblem *myProblem, double *U, double *F)
{
    int  size = myProblem->size, i;
    double h = 1.0/(size-1);
    double beta    = myProblem->beta;
    double epsilon = myProblem->epsilon;
    double zeta    = myProblem->zeta;

    F[0] = 0;
    for (i=1; i < size-1; i++)
        F[i] = (epsilon * (U[i-1] - 2 * U[i] + U[i+1])/(h*h) 
        - (1-zeta)*beta * (U[i+1] - U[i-1])/(2*h) 
        - zeta*beta*(U[i] - U[i-1])/(h));
	F[size-1] = 0;
}


# endif
# ifndef NOEULER

void advdiffUpdateEuler(femProblem *myProblem, double dt)
{
    int  size = myProblem->size, i;
    double *U = myProblem->U;
    double *F = malloc(sizeof(double) * size);

    advdiffComputeRightHandSide(myProblem,U,F);
    for (i=0; i < size; i++)
        U[i] += dt * F[i];

    free(F);
}

# endif
# ifndef NORK

void advdiffUpdateRungeKutta(femProblem *myProblem, double dt)
{
    int  size = myProblem->size, i,j;
    double *U = myProblem->U;
    double *Uold       = malloc(sizeof(double) * size);
    double *Upredictor = malloc(sizeof(double) * size);
    double *K = malloc(sizeof(double) * size);
    for (i=0; i < size; i++) {
        K[i] = 0.0;
        Uold[i] = U[i];  }

    const int    nStage   = 4;
    const double beta[4]  = {0.0,     0.5,   0.5, 1.0  };
    const double gamma[4] = {1.0/6, 2.0/6, 2.0/6, 1.0/6};

    for (j = 0; j < nStage; j++) {
        for (i=0; i < size; i++)
        	Upredictor[i] = Uold[i] + dt * beta[j] * K[i];
    	advdiffComputeRightHandSide(myProblem,Upredictor,K);
    	for (i=0; i < size; i++)
        	U[i] += dt * gamma[j] * K[i]; }

    free(Uold);
    free(Upredictor);
    free(K);
}

# endif