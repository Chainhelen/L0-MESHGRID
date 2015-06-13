#ifndef RECOV_H
#define RECOV_H

#include "glm.h"
#include "SharedStructures.h"
#include <gsl/gsl_matrix.h>
#include <map>
#include <vector>
#include <cmath>

class ReCov{
    public:
        IndexList **faceNeighborFace;
        IndexList **rFaceNeighborFace;
        IndexList **dFaceNeighborFace;
        IndexList **verticestindices;
        IndexList **verticesvindices;

        double *featureVector;
        double **faceNormal;
        double *vectorOfCoV;
        double *delta;
        double **weightCov;
        int length;
        GLMmodel *meshmodel;

    public:
        ReCov(GLMmodel *pmeshmodel,IndexList **pverticesvindices,IndexList **pverticestindices);
        ~ReCov();
        //面片r环领域
        void getFaceNeighborFace();
        void delFaceNeighborFace();
        int checkCommonVerticeOfNeighborFace(int afaceindex,int bfaceindex);
        void initDFaceNeighborFace();
        void delDFaceNeighborFace();
        void getFaceInsideFaces(IndexList **Head,IndexList **faceNeighborFace,int r);

        //面法向操作
        void initFaceNormal();
        void updateFaceNormal();
        void getFaceVectorByFaceIndex(double *vec, int findex);
        void delFaceNormal();

        //获得featurevector
        void getFeatureVector();
        void delFeatureVector();

        //计算D换领域协方差矩阵均值
        void computeMeanOfFeature(IndexList **dFaceNeighborFace,double *mean,int vindex);
        void computeCov(IndexList **dFaceNeighborFace,gsl_matrix *Cov,double *mean, int vindex);

        void getVectorOfCoV();
        void delVectorOfCov();

        //得到delta
        void getDelta();
        void delDelta();

        //得到权重
        void getWeightCov();
        void delWeightCov();

        void cholesky(gsl_matrix *A);
        void insertIndexList(IndexList **H,int place,std::map<int,char> &mp, int num);
        void insertIndexList(IndexList **H,int place,std::vector<int> &v, int num);

        void delIndexList(IndexList **Head,int num);


        void doCov();
};
#endif
