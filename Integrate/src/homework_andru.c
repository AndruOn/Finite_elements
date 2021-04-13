#include <stdio.h>
#include <math.h>
#include "glfem.h"


double integrate(double x[3], double y[3], double (*f) (double, double)){
    double xhi[3] = {1.0/6 , 1.0/6 , 2.0/3};
    double eta[3] = {1.0/6 , 2.0/3 , 1.0/6};

    double xLoc[3];
    double yLoc[3];

    for (size_t j = 0; j < 3; j++){
        xLoc[j] = x[0] * (1-xhi[j] - eta[j]) + x[1] * xhi[j] + x[2] * eta[j];
        yLoc[j] = y[0] * (1-xhi[j] - eta[j]) + y[1] * xhi[j] + y[2] * eta[j];;
    }
    

    glfemSetColor(GLFEM_BLACK); glfemDrawElement(x,y,3);
    glfemSetColor(GLFEM_BLUE);  glfemDrawNodes(x,y,3);
    glfemSetColor(GLFEM_RED);   glfemDrawNodes(xLoc,yLoc,3);
    
    //need to integrate
    double jacobien = fabs( 
        ( (x[1]-x[0])*(y[2]-y[0]) ) 
        - 
        ( (x[2]-x[0])*(y[1]-y[0]) ) 
    );
    /*
    double B[3] = {0.5 * jacob, 1.0/6 * jacob, 1.0/6 * jacob};
    double A[3][3] = {
        {1,1,1},
        {xLoc[0],xLoc[1],xLoc[2]},
        {yLoc[0],yLoc[1],yLoc[2]}
    };
    //TO FINNNND
    double W[3]= {1.0/6,1.0/6,1.0/6};//solve(A,B); //need to find slove systemvfun
    //TO FIINNNND
    */
   double W[3]= {1.0/6,1.0/6,1.0/6};
    double I = 0;
    for (size_t i = 0; i < 3; i++){
        I+= W[i] * f(xLoc[i],yLoc[i]);
    }

    return jacobien * I;
}



double integrateRecursive(double x[3], double y[3], double (*f)(double,double), int n)
{
    if(n == 0){
        return integrate(x,y,f);
    }
    double half_x[3] = {(x[2]+x[0])/2, (x[0]+x[1])/2, (x[1]+x[2])/2};
    double half_y[3] = {(y[2]+y[0])/2, (y[0]+y[1])/2, (y[1]+y[2])/2};
    //for (size_t i = 0; i < 3; i++) printf("halfx:%f halfy:%f\n",half_x[i],half_y[i]);

    double I = 0;
    //printf("///////////////:n:%d//////////////////\n",n);
    for (size_t i = 0; i < 3; i++){
        //printf("i:%d\n",i);
        double new_x[3] = {x[i], half_x[i], half_x[(i+1)%3]};
        double new_y[3] = {y[i], half_y[i], half_y[(i+1)%3]};
        //for(size_t j = 0; j < 3; j++) printf("x:%f y:%f\n",new_x[j],new_y[j]);
        //integrate(x,y,f);
        I += integrateRecursive(new_x,new_y,f,n-1);
    }
    //printf("i:3\n");
    double new_x[3] = {half_x[0], half_x[1], half_x[2]};
    double new_y[3] = {half_y[0], half_y[1], half_y[2]};
    //for(size_t j = 0; j < 3; j++) printf("x:%f y:%f\n",new_x[j],new_y[j]);
    //integrate(x,y,f);
    I+= integrateRecursive(new_x,new_y,f,n-1);
    return I;
}
