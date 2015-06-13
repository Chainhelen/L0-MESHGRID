#ifndef L0BYFACEVECTOR_H
#define L0BYFACEVECTOR_H
#include "L0.H"
#include "RECOVERYBYL2.H"

class L0ByFaceVector : public L0{
    public:
        Info **info;
        List **faceNeighborFace;
        List **faceSpreadNeighborFace;
        double **p;
        double **v;
        double **weightCov;
    public:
        L0ByFaceVector(GLMmodel *originmeshmodel,IndexList **originverticesvindices,IndexList **originverticestindices);
        ~L0ByFaceVector();


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
