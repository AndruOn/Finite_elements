#include <stdio.h>
#include <math.h>
#include "glfem.h"


double integrate(double x[3], double y[3], double (*f) (double, double)){
    double E[3] = {1.0/6 , 1.0/6 , 2.0/3};
    double N[3] = {1.0/6 , 2.0/3 , 1.0/6};

    
    double I = 0;
    double xLoc[3];
    double yLoc[3];

    for (size_t i = 0; i < 3; i++)
    {
        xLoc[i] = x[0]*(1-E[i]-N[i]) + x[1]*(E[i]) + x[2]*N[i];
        yLoc[i] = y[0]*(1-E[i]-N[i]) + y[1]*(E[i]) + y[2]*N[i];
    }
    

    
    glfemSetColor(GLFEM_BLACK); glfemDrawElement(x,y,3);
    glfemSetColor(GLFEM_BLUE);  glfemDrawNodes(x,y,3);
    glfemSetColor(GLFEM_RED);   glfemDrawNodes(xLoc,yLoc,3);

    double J = (x[1]-x[0])*(y[2]-y[0])-(x[2]-x[0])*(y[1]-y[0]);
    
    for (size_t i = 0; i < 3; i++)
    {
        I += 1.0/6*f(xLoc[i], yLoc[i]);
    }
    
    return fabs(J)*I;
}



double integrateRecursive(double x[3], double y[3], double (*f)(double,double), int n)
{

    double I =0;
    if(n == 0){
        return integrate(x,y,f);
    }


    else
    {
        double middlex[3];
        double middley[3];
        


        middlex[0] = (x[0] + x[1])/2;
        middlex[1] = (x[0] + x[2])/2;
        middlex[2] = (x[1] + x[2])/2;
        middley[0] = (y[0] + y[1])/2;
        middley[1] = (y[0] + y[2])/2;
        middley[2] = (y[1] + y[2])/2;

        double ptr1x[3] = {x[0], middlex[0], middlex[1]}; 
        double ptr2x[3] = {middlex[2], middlex[0], x[1]};
        double ptr3x[3] = {middlex[2], middlex[1], x[2]};
        double ptr1y[3] = {y[0], middley[0], middley[1]};
        double ptr2y[3] = {middley[2], middley[0], y[1]};
        double ptr3y[3] = {middley[2], middley[1], y[2]};

        
        
        I= integrateRecursive(ptr1x, ptr1y, f, n-1)  
        + integrateRecursive(ptr2x, ptr2y, f, n-1) 
        + integrateRecursive(ptr3x, ptr3y, f, n-1) 
        + integrateRecursive(middlex, middley, f, n-1); 
    }
    return I;



    
}
