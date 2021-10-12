#ifndef SIMPSONHEADERDEF
#define SIMPSONHEADERDEF

#include "AbstractQuadratureRule.hpp"

class Simpson:
    public AbstractQuadratureRule
{

    public:

        //Specialised constructor
        Simpson(double (*pFunction) (double), const double xmin, const double xmax,
                const int NoValues);

        double IntegrateFunction();

        double IntegrateRHSProduct
            (const int i, int npoints, Vector* pPoints);

        double IntegrateMatrixProduct
            (const int i, const int j, int npoints, Vector* pPoints);

    private:

        //Hidden default constructor
        Simpson();


};

#endif
