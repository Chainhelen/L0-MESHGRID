#ifndef L0BYFACENEIGHBOREDGE_H
#define L0BYFACENEIGHBOREDGE_H
#include "L0.H"

class L0ByFaceNeighborEdge : public L0{
    public:
        Info **info;
        List **faceNeighborFace;
        double *p;
        double *v;
        int *infoToFaceIndex;
    public:
        L0ByFaceNeighborEdge(GLMmodel *originmeshmodel,IndexList **originverticesvindices,IndexList **originverticestindices);
        ~L0ByFaceNeighborEdge();
        GLMmodel* doL0(double parpha, double pbeta, double plambda, int pmaxtimes);

        void chooseTypeOfNeighbor();
        void delFaceNeighborFace();
        void initInfo();
        void updateInfo();
        void delInfo();
};
#endif
