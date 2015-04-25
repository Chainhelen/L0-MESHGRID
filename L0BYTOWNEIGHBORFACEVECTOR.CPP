#include "L0BYTOWNEIGHBORFACEVECTOR.H"

//#define FIRSTFACENEIGHBORFACE
//#define SECONDFACENEIGHBORFACE
#define ALLFACENEIGHBORFACE


L0ByTowNeighborFaceVector::L0ByTowNeighborFaceVector(GLMmodel *originmeshmodel, IndexList **originverticesvindices, IndexList **originverticestindices) \
        :L0(originmeshmodel, originverticesvindices, originverticestindices)
{
    info = NULL;
    faceNeighborFace = NULL;
    faceSpreadNeighborFace = NULL;
    p = NULL;
    v = NULL;
}
L0ByTowNeighborFaceVector::~L0ByTowNeighborFaceVector()
{
    if(NULL != faceVector)
    {
        delVector(faceVector);
        faceVector = NULL;
    }
    if(v){
        for(int i = 0;i < 3;i++){
            delete []v[i];
            v[i] = NULL;
        }
        delete []v;
        v = NULL;
    }
    if(p){
        for(int i = 0;i < 3;i++){
            delete []p[i];
            p[i] = NULL;
        }
        delete []p;
        p = NULL;
    }
    delInfo();
    delFaceNeighborFace();
    delFaceSpreadNeighborFace();

    faceNeighborFace = NULL;
    faceSpreadNeighborFace = NULL;
}

GLMmodel* L0ByTowNeighborFaceVector::doL0(double parpha, double pbeta, double plambda, int pmaxtimes)
{
    arpha = parpha;
    beta = pbeta;
    lambda = plambda;
    maxtimes = pmaxtimes;

    return meshmodel;
}
void L0ByTowNeighborFaceVector::chooseTypeOfNeighbor()
{
#ifdef FIRSTFACENEIGHBORFACE
   initFaceFirstNeighborFace();
   faceNeighborFace = faceFirstNeighborFace;
#endif

#ifdef SECONDFACENEIGHBORFACE
   initFaceSecondNeighborFace();
   faceNeighborFace = faceSecondNeighborFace;
#endif

#ifdef ALLFACENEIGHBORFACE
   initFaceAllNeighborFace();
   faceNeighborFace = faceAllNeighborFace;
#endif
   initFaceSpreadNeighborFace();
}

void L0ByTowNeighborFaceVector::initInfo()
{

}
void L0ByTowNeighborFaceVector::updateInfo()
{

}
void L0ByTowNeighborFaceVector::delInfo()
{

}
void L0ByTowNeighborFaceVector::delFaceNeighborFace()
{

}
void L0ByTowNeighborFaceVector::initFaceSpreadNeighborFace()
{

}
void L0ByTowNeighborFaceVector::delFaceSpreadNeighborFace()
{

}
void L0ByTowNeighborFaceVector::recoveryVerticesByFaceVector()
{

}
