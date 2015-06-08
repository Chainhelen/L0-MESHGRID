#include "CovVerNormal.h"
#include "RecoveryByVerNormalL0.h"

CovVerNormal::CovVerNormal(GLMmodel *originmeshmodel, IndexList **originverticesvindices, IndexList **originverticestindices) \
        :L0(originmeshmodel, originverticesvindices, originverticestindices)
{
	verWeight = NULL;
}

CovVerNormal::~CovVerNormal()
{
    delVerWeight();
}

void CovVerNormal::initVerWeight()
{
    int i;
    verWeight = new double*[3];

    for(i = 0;i < 3;i++){
        verWeight[i] = new double[(int)meshmodel->numvertices];
    }
}

void CovVerNormal::delVerWeight()
{
    int i;
    if(NULL != verWeight){
        for(i = 0;i < 3;i++){
            delete[] verWeight[i];
            verWeight[i] = NULL;
        }
        delete[] verWeight;
        verWeight = NULL;
    }
}

void CovVerNormal::recoveryVerticesByVerNormal()
{
    //printf("vernum = %d\ntrinum = %d\n",(int)meshmodel->numvertices,(int)meshmodel->numtriangles);
    initVerNeighborVer();
    RecoveryByVerNormalL0 myrecoverybyl0(meshmodel, verNeighborVer,verVector);
    delVerNeighborVer();

    getEdge();
    initVerSpreadNeighborVer();

    myrecoverybyl0.getParameter(verSpreadNeighborVer, edgeNum, edge);
    myrecoverybyl0.slove();

    delVerSpreadNeighborVer();
    delEdge();
}
