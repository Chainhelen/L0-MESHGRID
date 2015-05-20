#ifndef RECOVERYBYVERNORMALL0_H
#define RECOVERYBYVERNORMALL0_H
#include "L0.H"

class RecoveryByVerNormalL0{
    public:
        GLMmodel *meshmodel;
        IndexList **verticesvindices,**verticestindices;
        double **vervector;

        Info **info;
        Edge *edge;
        double *p;
        double *v;

        List **verRelation;
        List **relation;
        int edgeCnt;
        int arinfocnt;

    public:
        RecoveryByVerNormalL0(GLMmodel *pmeshmodel,IndexList **pverticesvindices,IndexList **pverticestindices,double **pvervector);
		~RecoveryByVerNormalL0();

        void getParameter(List **pVerRelation,int pEdgeCnt, Edge *pEdge);
        void initRelation();
        void initInfo();
        void slove();
};
#endif
