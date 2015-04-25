#ifndef L0_H
#define L0_H

#include "glm.h"
#include <cmath>
#include <stdio.h>
#include <algorithm>
#include <vector>
#include <iterator>
#include "SharedStructures.h"
#include "SubSolving.h"

using namespace std;


class L0{
    public:
        GLMmodel *meshmodel;
        IndexList **verticesvindices;
        IndexList **verticestindices;
        double arpha, beta, lambda;

        double *model;
        double **faceVector;
        double **verVector;
        int maxtimes;
        int edgeNum;
        int arinfocnt, infocnt;

        List **verNeighborVer;
        List **verSpreadNeighborVer;
        List **faceFirstNeighborFace;
        List **faceSecondNeighborFace;
        List **faceAllNeighborFace;
        List **faceNeighborEdge;

        Edge *edge;

    public:
        L0(GLMmodel *originmeshmodel,IndexList **originverticesvindices,IndexList **originverticestindices);
        ~L0();

        void initModel();
        //点的领域包括自己
        //点的一环领域
        void initVerNeighborVer();
        void delVerNeighborVer();

        //点的二环领域
        void initVerSpreadNeighborVer();
        void delVerSpreadNeighborVer();

        //面的领域包括自己
        //面的一类领域点
        void initFaceFirstNeighborFace();
        void delFaceFirstNeighborFace();

        //面的二类领域点
        void initFaceSecondNeighborFace();
        void delFaceSecondNeighborFace();

        //面的所有领域点
        void initFaceAllNeighborFace();
        void delFaceAllNeighborFace();

        //面领域与周围面片的边内积
        void initFaceNeighborEdge(List **faceNeighborFace);
        void delFaceNeighborEdge();

        //init法相不具有初始法相
        void initFaceVector();
        void initVerVector();
        void updateFaceVector();
        void updateVerVectorByArea();

        //删除点，面法相
        void delVector(double **pvector);


//判断类函数
        int checkCommonEdge(int x, int y, int commonindex);
        int checkCommonVerticeOfNeighborFace(int afaceindex, int bfaceindex);

//辅助类函数
        double getPointDistance(int x, int y);
        double getPointAngle(double z, double x, double y);
        double getPointFactor(int index1, int index2, int index3, int index4);
        double getAreaByIndex(int findex);
        double getDistanceFaceCenterToFaceCenter(int afindex,int bfindex);
        double getAngleByFaceIndex(int afindex,int bfindex);
        double getSignalByFaceIndex(int afindex,int bfindex);
	    void getFaceVectorByFaceIndex(double *vec, int findex);

        void getEdge();
        void delEdge();

        Edge *createEdgeNode();
        List *createListNode();
        void insertList(List **relation, int place, int *a, int num);
        void insertList(List **relation, int place, std::vector<int> v, int num);
};

#endif L0_H
