#include"fem.h"



#ifndef NOASSEMBLE


void femIterativeSolverAssemble(femIterativeSolver* mySolver, double *Aloc, double *Bloc, double *Uloc, int *map, int nLoc)
{
    int iLoc,j;
    int myRow ;
    if (mySolver->iter==0)
    {
        for (iLoc = 0; iLoc < nLoc; iLoc++) 
        {
            myRow = map[iLoc];
            mySolver->R[myRow] += Bloc[iLoc];
            for(j = 0; j < nLoc; j++) 
            {
                mySolver->R[myRow] -= Aloc[iLoc*nLoc+j]*Uloc[j]; 
            }
        }
    }
    else
    {   
        for (iLoc = 0; iLoc < nLoc; iLoc++) 
        {
            for(j = 0; j < nLoc; j++) 
            {
                myRow = map[iLoc];
                mySolver->S[myRow]+=Aloc[iLoc*nLoc+j]*mySolver->D[map[j]];//Aloc[iLoc*nLoc+j]*Uloc[j];
            }
        }
    }
               
}
/*{
    int iLoc,j;
    for (iLoc = 0; iLoc < nLoc; iLoc++) { 
        int myRow = map[iLoc];
        mySolver->R[myRow] += Bloc[iLoc];
        
        for(j = 0; j < nLoc; j++) {
            mySolver->R[myRow] -= Aloc[iLoc*nLoc+j]*Uloc[j];
        }
    }
}
*/

#endif
#ifndef NOITERATE


double *femIterativeSolverEliminate(femIterativeSolver *mySolver)
{   mySolver->iter++;
    double error = 0.0; int i;
    double alpha = 0.2;
    double beta=0;
    if (mySolver->iter==1)
    {
        for (i=0; i < mySolver->size; i++)
        {
            mySolver->D[i]=mySolver->R[i]; 
            mySolver->X[i] = 0.0; 
        }
    }
    else
    {
        //alpha
        alpha=0;
        double inter1=0;
        double inter2=0;
        for (i=0; i < mySolver->size; i++) 
        {
            inter1+=mySolver->R[i]*mySolver->R[i];
            inter2+=mySolver->S[i]*mySolver->R[i];
        }
        alpha=inter1/inter2;
        inter2=0;
        //incrémente R
        for (i=0; i < mySolver->size; i++) 
        {   
            mySolver->R[i]-=alpha*mySolver->S[i];
        }
        //Calcul de Beta
        for (i=0; i < mySolver->size; i++) 
        {
            inter2+=mySolver->R[i]*mySolver->R[i];
            mySolver->X[i] = mySolver->D[i]*alpha; 
        }
        beta=inter2/inter1;
        //recacul de D
        for (i=0; i < mySolver->size; i++) 
        {   
            mySolver->D[i]=mySolver->R[i]+beta*mySolver->D[i];
        }
        //relttre a 0
        for (i=0; i < mySolver->size; i++) 
        {
            mySolver->S[i]=0;
        }
        mySolver->error=sqrt(inter1);
    }
        return(mySolver->X);

/*
    mySolver->iter++;
    double error = 0.0; int i;
    double alpha = 0.2;
    for (i=0; i < mySolver->size; i++) {
        error += (mySolver->R[i])*(mySolver->R[i]);
        mySolver->X[i] = mySolver->R[i]*alpha; 
        mySolver->R[i] = 0.0; }
        
    mySolver->error = sqrt(error);
    return(mySolver->X);
    */
    
}

#endif