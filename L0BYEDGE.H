#ifndef L0BYEDGE_H
#define L0BYEDGE_H
#include "L0.H"

class L0ByEdge : public L0{
    public:
        Info **info;
        double **p;
        double **v;
    public:
        L0ByEdge(GLMmodel *originmeshmodel,IndexList **originverticesvindices,IndexList **originverticestindices);
        ~L0ByEdge();
        GLMmodel* doL0(double parpha, double pbeta, double plambda, int pmaxtimes);

        void initInfo();
        void updateInfo();
        void delInfo();
        void autoGetParameter();
};
#endif
