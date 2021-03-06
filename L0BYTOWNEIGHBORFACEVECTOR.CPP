#include "L0BYTOWNEIGHBORFACEVECTOR.H"

//#define FIRSTFACENEIGHBORFACE
#define SECONDFACENEIGHBORFACE
//#define ALLFACENEIGHBORFACE


L0ByTowNeighborFaceVector::L0ByTowNeighborFaceVector(GLMmodel *originmeshmodel, IndexList **originverticesvindices, IndexList **originverticestindices) \
        :L0(originmeshmodel, originverticesvindices, originverticestindices)
{
    info = NULL;
    faceNeighborFace = NULL;
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

    faceNeighborFace = NULL;
}

GLMmodel* L0ByTowNeighborFaceVector::doL0(double parpha, double pbeta, double plambda, int pmaxtimes)
{
    int i,j,k;
    arpha = parpha;
    beta = pbeta;
    lambda = plambda;
    maxtimes = pmaxtimes;
    double tsum[3],sum;
    int cc = 1;

    initFaceVector();
    updateFaceVector();

    chooseTypeOfNeighbor();
    initInfo();
    updateInfo();

    getPV();


    SubSolving s_l0(faceNeighborFace, info, infocnt, arinfocnt, (int)meshmodel->numtriangles);
    s_l0.init();

    while(cc <= pmaxtimes ){
        updateInfo();
        for(i = arinfocnt;i < infocnt;i++){
            for(j = 0;j < 3;j++){
                tsum[j] = 0.0;
            }
            sum = 0.0;
            for(j = 0;j < info[i][0].cnt;j++){
                for(k = 0;k < 3;k++){
                    tsum[k] += info[i][j].w * faceVector[k][info[i][j].data];
                }
            }
            for(j = 0;j < 3;j++){
                sum += tsum[j] * tsum[j];
            }
            if(sum <= lambda / beta){
                for(j = 0;j < 3;j++){
                    p[j][i + (int)meshmodel->numtriangles] = 0.0;
                }
            }else{
                for(j = 0;j < 3;j++){
                    p[j][i + (int)meshmodel->numtriangles] = beta * tsum[j] ;
                }
            }
        }
        //传入beta,arpha
        s_l0.getParameter(p[0], v[0], info, beta, arpha);
        s_l0.update();

        for(i = 0;i < 3;i++){
            //            传入p,v

            s_l0.getParameter(p[i], v[i], info, beta, arpha);
            s_l0.slove();
            for(j = 0;j < (int)meshmodel->numtriangles;j++){
                faceVector[i][j] = v[i][j];
            }
        }

        for(i = 0;i < (int)meshmodel->numtriangles;i++){
            sum = 0.0;
            for(j = 0;j < 3;j++){
                sum += faceVector[j][i] * faceVector[j][i];
            }
            sum = sqrt(sum);
            sum = sum > 1e-12 ? sum : 1e-12;
            for(j = 0;j < 3;j++){
                faceVector[j][i] = faceVector[j][i] / sum;
            }

        }

        cc++;
        beta = beta * sqrt(2);
        arpha = arpha / 2;
    }
    recoveryVerticesByFaceVector();

    delVector(faceVector);
    faceVector = NULL;

    for(i = 0;i < (int)meshmodel->numvertices;i++){
        for(j = 0;j < 3;j++){
            meshmodel->vertices[3 * (i + 1) + j] = model[3 * i + j];
        }
    }

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
}

void L0ByTowNeighborFaceVector::initInfo()
{

    int i, j;

    infocnt = 0;
    arinfocnt = 0;
    List *faceneighborfacetail = NULL;
    vector<int>v;

    for(i  = 0;i < (int)meshmodel->numtriangles;i++){
        faceneighborfacetail = faceNeighborFace[i];
        while(faceneighborfacetail){
            if(i < faceneighborfacetail->data){
                v.push_back(i);
                v.push_back(faceneighborfacetail->data);
            }
            faceneighborfacetail = faceneighborfacetail->next;
        }
    }

    infocnt = (int)v.size() / 2;

    info = new Info*[(int)infocnt];
    if(NULL == info){
        printf("error at L0ByTowNeighborFaceVector::initInfo  in info\n");
        return ;
    }

    for(i = 0;i < infocnt;i++){
        info[i] = new Info[2];
        if(NULL == info[i]){
            printf("error at L0ByTowNeighborFaceVector::initInfo in info[%d]\n",i);
            return ;
        }
        for(j  = 0;j < 2;j++){
            info[i][j].cnt = 2;
            info[i][j].data = v[2 * i + j];
            info[i][j].w = 0.0;
        }
    }
}

void L0ByTowNeighborFaceVector::updateInfo()
{
    int i;
    for(i = arinfocnt;i < infocnt;i++){
        info[i][0].w = 1.0;
        info[i][1].w = -1.0;
    }
}
void L0ByTowNeighborFaceVector::delInfo()
{
    int i;
    if(info){
        for(i = 0;i < infocnt;i++){
            delete  []info[i];
            info[i] = NULL;
        }
        delete []info;
        info = NULL;
    }
}

void L0ByTowNeighborFaceVector::getPV()
{
    int i, j;
    p = new double*[3];
    v = new double*[3];

    if(NULL == p || NULL == v){
        printf("error in L0ByTowNeighborFaceVector::doL0 p v\n");
        return ;
    }

    for(i = 0;i < 3;i++){
        p[i] = new double[(int)meshmodel->numtriangles + infocnt];
        v[i] = new double[(int)meshmodel->numtriangles];

        if(NULL == p[i] || NULL == v[i]){
            printf("error in L0ByTowNeighborFaceVector::doL0 p[%d] v[%d]\n",i,i);
            return ;
        }
    }
    for(i = 0;i < 3;i++){
        for(j = 0;j < (int)meshmodel->numtriangles;j++){
            p[i][j] = faceVector[i][j];
            v[i][j] = 0;
        }
    }
    for(i = 0;i < 3;i++){
        for(j = (int)meshmodel->numtriangles;j < (int)meshmodel->numtriangles + infocnt;j++){
            p[i][j] = 0;
        }
    }
}

void L0ByTowNeighborFaceVector::delFaceNeighborFace()
{
#ifdef FIRSTFACENEIGHBORFACE
   delFaceFirstNeighborFace();
#endif

#ifdef SECONDFACENEIGHBORFACE
   delFaceSecondNeighborFace();
#endif

#ifdef ALLFACENEIGHBORFACE
   delFaceAllNeighborFace();
#endif
   faceNeighborFace = NULL;
}

void L0ByTowNeighborFaceVector::recoveryVerticesByFaceVector()
{
    RecoveryByL2 myRecovery(meshmodel, verticesvindices, verticestindices);
    myRecovery.changeModel(model);
    myRecovery.changeFaceVector(faceVector);
    myRecovery.doL2(2);
}
