#ifndef L0BYTOWNEIGHBORFACEVECTOR_H
#define  L0BYTOWNEIGHBORFACEVECTOR_H
#include "L0.H"
#include "RECOVERYBYL2.H"

class L0ByTowNeighborFaceVector : public L0{
    public:
        Info **info;
        List **faceNeighborFace;
        double **p;
        double **v;
        double **weightCov;
    public:
        L0ByTowNeighborFaceVector(GLMmodel *originmeshmodel,IndexList **originverticesvindices,IndexList **originverticestindices);
        ~L0ByTowNeighborFaceVector();


        void useWeightCov(double **pWeightCov);
        GLMmodel* doL0(double parpha, double pbeta, double plambda, int pmaxtimes);

        void chooseTypeOfNeighbor();

        void initInfo();
        void updateInfo();
        void delInfo();
        void getPV();

        void delFaceNeighborFace();
        void recoveryVerticesByFaceVector();
};
#endif
