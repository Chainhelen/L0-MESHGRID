
//Cov.h

#include "glm.h"
#include "SharedStructures.h"
#include <gsl/gsl_matrix.h>

#pragma once

class Cov
{
public:

	int length;                                                 //特征向量维数
	double *featureVector;                                      //特征向量
	IndexList **verOfNeighborVer;                               //每个顶点相邻点的索引，第一个节点存放邻点的个数
	IndexList **trianglesOfNeighborVer;                         //每个顶点相邻三角形的索引 ，第一个节点存放相邻三角形的个数
	double *facetnorm, *vernormal, *laplaciancoordinates;

	int n;                                                      //顶点个数；
	//权值
	double weightd, weightn, weightl;

	double *vectorOfCoV;

	double **weightCov, *delta;

	Cov();
	~Cov();

	void FindNeighborVerTriangles(GLMmodel *meshmodel); 
	void ComputeFaceNormal(GLMmodel *meshmodel);
	void ComputeVertexNormal(GLMmodel *meshmodel);

	//拉普拉斯坐标
	bool IsInTriangle(GLMmodel *pmesh, int indexv1, int indexv2,
		int *indexv3, int indext);
	void ComputeCotangentWeight(double *pcotangentweight, int indexv, int *indexcw, GLMmodel *pmesh, 
		IndexList **pverOfNeighborVer, IndexList **ptrianglesOfNeighborVer);
	void ConstructLaplacianCoordinatesCotangentWeight(double *plc, double *pcw, int *indexcw, int indexv,
		IndexList **pverOfNeighborVer, float *vertices);
	void ComputeLaplaciancoordinates(GLMmodel *pmesh, IndexList **pverOfNeighborVer, IndexList **ptrianglesOfNeighborVer);

	void GetFeatureVector(GLMmodel *meshmodel);
	void ReturnNumofVerofMeshodel(GLMmodel *meshmodel);

	//计算n环邻域
	void AddVertex(GLMmodel *meshmodel, IndexList **L);
	void ComputeNrneighbor(GLMmodel *meshmodel, IndexList **H, int r);

	void ComputeMeanOfFeature(IndexList **neighbor, double *mean, int indexv);
	void ComputeCov(IndexList **neighbor, gsl_matrix *cov, double *mean, int indexv);
	void Cholesky(gsl_matrix * A);
	void GetVectorOfCoV(GLMmodel *meshmodel);

	void Computedelta(GLMmodel *meshmodel);
	void ComputeWeightCov(GLMmodel *meshmodel);
	void Denoising(GLMmodel *meshmodel);

	void Deletes(IndexList **pverOfNeighborVer, IndexList **ptrianglesOfNeighborVer, double **pweightCov, int pn);

};

