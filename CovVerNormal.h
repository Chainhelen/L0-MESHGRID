#ifndef COVVERNORMAL_H
#define COVVERNORMAL_H
#include "L0.H"

class CovVerNormal : public L0{
    public:
		double **verWeight;
    public:
        CovVerNormal(GLMmodel *originmeshmodel,IndexList **originverticesvindices,IndexList **originverticestindices);
        ~CovVerNormal();
        void initVerWeight();
        void delVerWeight();
        void recoveryVerticesByVerNormal();
};
#endif
