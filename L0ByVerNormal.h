#ifndef L0BYVERNORMAL_H
#define L0BYVERNORMAL_H
#include "L0.H"

class L0ByVerNormal : public L0{
    public:
        Info **info;
        double **p;
        double **v;
        GLMmodel *guidedmesh;
        double **guidedmeshVerNormal;
        double **guidedmeshFaceNormal;

    public:
        L0ByVerNormal(GLMmodel *originmeshmodel,IndexList **originverticesvindices,IndexList **originverticestindices);
        ~L0ByVerNormal();

        void getGuidedmesh(GLMmodel *pGuidedmesh);
        void initGuidedmeshFaceNormal();
        void initGuidedmeshVerNormal();

        GLMmodel* doL0(double parpha, double pbeta, double plambda, int pmaxtimes);

        void getPV();
        void initInfo();
        void updateInfo();
        void delInfo();

        void recoveryVerticesByVerNormal();
};
#endif
