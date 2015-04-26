#include "StdAfx.h"
#include "stdlib.h"
#include "SharedStructures.h"
#include "SolveUnSymmetricSparseSystem.h"
#include "windows.h"
#include "string.h"
#include <algorithm>
#include "glm.h"

#ifndef SUBSOLVING_H
#define SUBSOLVING_H

struct List{
    List *next;
    int data;
};

struct Info{
    int data;
    int cnt;
    double w;
};

struct Edge{
    int p[4];
	int f[2];
    Edge *next;
};

template <class T>
void removeList(T *Head)
{
    if(NULL == Head)
        return;
    T *pNode = Head;
    T *pDel = NULL;

    while(pNode->next)
    {
        pDel = pNode->next;
        pNode->next = pDel->next;
        delete pDel;
    }
    delete pNode;
    pNode = NULL;
    pDel = NULL;
}
class SubSolving{
    public:

        List **relation;
        Info **info;

        SparseMatrix **AT;
        SparseMatrix **C;
        double *p;
        double *pb;
        double *v;
        double beta;
        double arpha;
        int infocnt;
        int arinfocnt;
        int row;
        MKL mkl;

    public:
        SubSolving(List **prelation,Info **pinfo,int pinfocnt,int parinfcnt, int prow);
        ~SubSolving();
        void init();
        void update();
        void slove();
        void getATMem();
        void getCMem();
        void getParameter(double *pp,double *pv,Info **pinfo,double pbeta,double parpha);
        void removeSparseMatrix(SparseMatrix *Head);

        void updateC();
        void updateAT();
        double Sparsesparse(SparseMatrix *A,SparseMatrix *B);
};

#endif SUBSOLVING_H
