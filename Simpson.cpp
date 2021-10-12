

#include "Simpson.hpp"


//Specialised constructor
Simpson::Simpson(double (*pFunction) (double), const double xmin, const double xmax,
const int NoValues)
{
    SetInterval(xmin, xmax);
    mpFunction = pFunction;
    mNoValues = NoValues;
}


double Simpson::IntegrateFunction()
{

    double stepsize = (mXmax-mXmin)/ (mNoValues-1);

    //Output variable
    double IntValue= 0.0;

    for (int k=0; k<mNoValues;k++)
    {
        if (k == 0) // creates weights depending on the value of k
        {
            IntValue += ((mpFunction)(mXmin + k*stepsize))*(stepsize/3.0);
        }
        else if (k % 2 == 1) // Tests if k is odd
        {
            IntValue += ((mpFunction)(mXmin + k*stepsize))*(4.0*stepsize/3.0);
        }
        else
        {
            IntValue += ((mpFunction)(mXmin + k*stepsize))*(2.0*stepsize/3.0);
        }

    }


    return IntValue;
}

double Simpson::IntegrateRHSProduct
(const int i, int npoints, Vector* pPoints)
{

    double stepsize = (mXmax-mXmin)/ (mNoValues-1);


    //Output variable
    double IntValue = 0.0;

    for (int k=0; k<mNoValues;k++)
    {
        if (k == 0) // creates weights depending on the value of k
        {
            //Now also multiplies by the Lagrange polynomial at the step point
            IntValue += ((mpFunction)(mXmin + k*stepsize))*(stepsize/3.0)
            *EvaluateLagrangeBasis(mXmin + k*stepsize,i,npoints,pPoints);
        }
        else if (k % 2 == 1) // Tests if k is odd
        {
            IntValue += ((mpFunction)(mXmin + k*stepsize))*(4.0*stepsize/3.0)
            *EvaluateLagrangeBasis(mXmin + k*stepsize,i,npoints,pPoints);
        }
        else
        {
            IntValue += ((mpFunction)(mXmin + k*stepsize))*(2.0*stepsize/3.0)
            *EvaluateLagrangeBasis(mXmin + k*stepsize,i,npoints,pPoints);
        }

    }


    return IntValue;
}


double Simpson::IntegrateMatrixProduct
(const int i, const int j, int npoints, Vector* pPoints)
{


    double stepsize = (mXmax-mXmin)/ (mNoValues-1);


    //Output variable
    double IntValue = 0.0;

    for (int k=0; k<mNoValues;k++)
    {
        if (k == 0) // creates weights depending on the value of k
        {
            //Now just multiplies Lagrange polynomials at i and j
            IntValue += EvaluateLagrangeBasis(mXmin + k*stepsize,j,npoints,pPoints)*(stepsize/3.0)
            *EvaluateLagrangeBasis(mXmin + k*stepsize,i,npoints,pPoints);
        }
        else if (k % 2 == 1) // Tests if k is odd
        {
            IntValue += EvaluateLagrangeBasis(mXmin + k*stepsize,j,npoints,pPoints)*(4.0*stepsize/3.0)
            *EvaluateLagrangeBasis(mXmin + k*stepsize,i,npoints,pPoints);
        }
        else
        {
            IntValue += EvaluateLagrangeBasis(mXmin + k*stepsize,j,npoints,pPoints)*(2.0*stepsize/3.0)
            *EvaluateLagrangeBasis(mXmin + k*stepsize,i,npoints,pPoints);
        }

    }


    return IntValue;
}
