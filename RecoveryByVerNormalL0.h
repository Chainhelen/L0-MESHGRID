#ifndef RECOVERYBYVERNORMALL0_H
#define RECOVERYBYVERNORMALL0_H
#include "L0.H"

class RecoveryByVerNormalL0{
    public:
        GLMmodel *meshmodel;
        List **verticesortvindices;
        double **vervector;

        Info **info;
        Edge *edge;
        double *p;
        double *v;
        int infocnt;
        int arinfocnt;

        List **verRelation;
        List **relation;
        int edgeCnt;

    public:
        RecoveryByVerNormalL0(GLMmodel *pmeshmodel,List **verticesortvindices,double **pvervector);
		~RecoveryByVerNormalL0();

        void getParameter(List **pVerRelation,int pEdgeCnt, Edge *pEdge);
        void initRelation();
        void initInfo();
        void slove();
};
#endif
