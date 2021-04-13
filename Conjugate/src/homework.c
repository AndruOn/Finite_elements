
#include"fem.h"



#ifndef NOASSEMBLE


void femIterativeSolverAssemble(femIterativeSolver* mySolver, double *Aloc, double *Bloc, double *Uloc, int *map, int nLoc)
{
    int i,j;
    if(mySolver->iter == 0){
        for (i = 0; i < nLoc; i++) { 
            int myRow = map[i];
            mySolver->R[myRow] += Bloc[i];
            for(j = 0; j < nLoc; j++) {
                mySolver->R[myRow] -= Aloc[i*nLoc+j]*Uloc[j]; 
            }
        }
    }else{
        double* D = mySolver->D;
        for (i = 0; i < nLoc; i++) { 
            int myRow = map[i];
            
            for(j = 0; j < nLoc; j++) {
                int myColumn = map[j];
                mySolver->S[myRow] += Aloc[i*nLoc+j]*D[myColumn]; 
            }
        }
    }
}

#endif
#ifndef NOITERATE


double *femIterativeSolverEliminate(femIterativeSolver *mySolver)
{
    mySolver->iter++;
    double error = 0.0; int i;
    double alpha;
    double beta = 0;
    int size = mySolver->size;

    //Premiere itération on setup tout et pose d_0 et r_0
    if(mySolver->iter == 1){
        for (i = 0; i < size; i++){
            mySolver->X[i] = 0.0; 
            mySolver->D[i] = mySolver->R[i];
        }
    }
    //iter>1 on applique l'algo de gradient conjugués
    else{
        double* R = mySolver->R;

        //Calculate alpha
        double r_ksquarred = 0.0; 
        double den_alpha = 0.0;
        for (i = 0; i < size; i++){
            r_ksquarred += R[i] * R[i];
            den_alpha += mySolver->S[i] * R[i];
        }
        alpha = r_ksquarred/den_alpha;
        
        //r^(k+1) et x^(k+1)
        for (i = 0; i < size; i++){
            R[i] -= alpha* mySolver->S[i];
            mySolver->X[i] = mySolver->D[i]*alpha; 
        }

        //Beta
        double r_kplus1_squarred = 0.0;
        for (i = 0; i < mySolver->size; i++){
            r_kplus1_squarred += R[i] * R[i];
        }
        beta = r_kplus1_squarred / r_ksquarred;

        //d^k+1
        for (i = 0; i < size; i++){
            mySolver->D[i] = R[i] + beta* mySolver->D[i];
        }

        //Réinitialise S a vecteur nul
        for (i = 0; i < size; i++){
            mySolver->S[i] = 0.0;
        }

        //Calcul erreur
        mySolver->error = sqrt(r_ksquarred);
        
    }
    return(mySolver->X);
}

#endif




