#ifndef L0BYFACENORMAL_H
#define L0BYFACENORMAL_H
#include "L0.H"
#include "RECOVERYBYL2.H"

class L0ByFaceNormal : public L0{
    public:
        Info **info;
        List **faceNeighborFace;
        List **faceSpreadNeighborFace;
        double **p;
        double **v;
        double **weightCov;
    public:
        L0ByFaceNormal(GLMmodel *originmeshmodel,IndexList **originverticesvindices,IndexList **originverticestindices);
        ~L0ByFaceNormal();


        void useWeightCov(double **pWeightCov);
        GLMmodel* doL0(double parpha, double pbeta, double plambda, int pmaxtimes);

        void chooseTypeOfNeighbor();
        void initInfo();
        void updateInfo();
        void delInfo();
        void delFaceNeighborFace();
        void initFaceSpreadNeighborFace();
        void delFaceSpreadNeighborFace();
        void recoveryVerticesByFaceVector();
};
#endif
