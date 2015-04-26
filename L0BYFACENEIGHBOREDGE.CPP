#include "L0BYFACENEIGHBOREDGE.H"

//#define FIRSTFACENEIGHBORFACE
//#define SECONDFACENEIGHBORFACE
#define ALLFACENEIGHBORFACE

//#define TESTZHENG_LOCAL
//#define TEST2
#define  MAX_TMP 1

L0ByFaceNeighborEdge::L0ByFaceNeighborEdge(GLMmodel *originmeshmodel, IndexList **originverticesvindices, IndexList **originverticestindices) \
        :L0(originmeshmodel, originverticesvindices, originverticestindices)
{
    info = NULL;
    faceNeighborFace = NULL;
    infoToFaceIndex = NULL;
    p = NULL;
    v = NULL;
}

L0ByFaceNeighborEdge::~L0ByFaceNeighborEdge()
{
    if(p){
        delete []p;
        p = NULL;
    }
    if(v){
        delete []v;
        v = NULL;
    }
    if(NULL ==  infoToFaceIndex){
        delete []infoToFaceIndex;
        infoToFaceIndex = NULL;
    }

    delInfo();
    delFaceNeighborFace();

    info = NULL;
    faceNeighborFace = NULL;
}

GLMmodel* L0ByFaceNeighborEdge::doL0(double parpha, double pbeta, double plambda, int pmaxtimes)
{
    arpha = parpha;
    beta = pbeta;
    lambda = plambda;
    maxtimes = pmaxtimes;
    int cc = 1;

    chooseTypeOfNeighbor();
    initInfo();

    while(cc <= pmaxtimes){

    }
    return NULL;
}

void L0ByFaceNeighborEdge::chooseTypeOfNeighbor()
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

void L0ByFaceNeighborEdge::delFaceNeighborFace()
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

void L0ByFaceNeighborEdge::initInfo()
{
    List *faceneighborfacetail = NULL;
    infocnt = 0;
    arinfocnt = 0;

    int a, b, c, i, j, k, fIndex;
    int num;
    vector<int>v;

    for(i  = 0 ;i < (int)meshmodel->numtriangles;i++){
        faceneighborfacetail = faceNeighborFace[i];
        while(faceneighborfacetail){
            infocnt++;
            faceneighborfacetail = faceneighborfacetail->next;
        }
        infocnt--;
    }
    infocnt *= 2;

    info = new Info*[infocnt];
    infoToFaceIndex = new int[infocnt];
    if(NULL == info){
        printf("error in L0BYFaceNeighborEdge\n");
        return ;
    }

    num = 0;
    for(i = 0; i < (int)meshmodel->numtriangles; i++){
        v.clear();
        faceneighborfacetail = faceNeighborFace[i];
        while(faceneighborfacetail){
            fIndex = faceneighborfacetail->data;
            if(i != fIndex){
                if(1 == checkCommonVerticeOfNeighborFace(i, fIndex)){
                    c = b = a = -1;
                    for(j = 0;j < 3 && a == -1;j++){
                        for(k = 0;k < 3 && a == -1;k++){
                            if(meshmodel->triangles[fIndex].vindices[j] == meshmodel->triangles[i].vindices[k]){
                                a = meshmodel->triangles[fIndex].vindices[j];
                            }
                        }
                    }
                    for(j = 0;j < 3;j++){
                        if(meshmodel->triangles[fIndex].vindices[j] != a){
                            if(b != -1){
                                b = meshmodel->triangles[fIndex].vindices[j];
                            }else{
                                c = meshmodel->triangles[fIndex].vindices[j];
                            }
                        }
                    }
                }
                if(2 == checkCommonVerticeOfNeighborFace(i, fIndex)){
                    c = b = a =-1;
                    for(j = 0;j < 3 && a == -1;j++){
                        for(k = 0;k < 3 && a == -1;k++){
                            if(meshmodel->triangles[fIndex].vindices[j] == meshmodel->triangles[i].vindices[k]){
                                break;
                            }
                        }
                        if(k == 3){
                            a = meshmodel->triangles[fIndex].vindices[j];
                        }
                    }
                    for(j = 0;j < 3;j++){
                        if(meshmodel->triangles[fIndex].vindices[j] != a){
                            if(b != -1){
                                b = meshmodel->triangles[fIndex].vindices[j];
                            }else{
                                c = meshmodel->triangles[fIndex].vindices[j];
                            }
                        }
                    }
                }
                info[num] = new Info[6];
                if(NULL == info[num]){
                    printf("error at L0ByFaceneighborFace--->info[%d]\n",num);
                    return ;
                }
                for(j = 0;j < 3;j++){
                    info[num][j].w = 0.0;
                    info[num][j].cnt = 6;
                    info[num][j].data = 3 * a + j;

                    info[num][3 + j].w = 0.0;
                    info[num][3 + j].cnt = 6;
                    info[num][3 + j].data = 3 * b + j;
                }
                infoToFaceIndex[num] = i;
                num++;

                info[num] = new Info[6];
                if(NULL == info[num]){
                    printf("error at L0ByFaceneighborFace--->info[%d]\n",num);
                    return ;
                }
                for(j = 0;j < 3;j++){
                    info[num][j].w = 0.0;
                    info[num][j].cnt = 6;
                    info[num][j].data = 3 * a + j;

                    info[num][3 + j].w = 0.0;
                    info[num][3 + j].cnt = 6;
                    info[num][3 + j].data = 3 * c + j;
                }
                infoToFaceIndex[num] = i;
                num++;
            }
            faceneighborfacetail = faceneighborfacetail->next;
        }
    }
    faceneighborfacetail = NULL;
}

void L0ByFaceNeighborEdge::updateInfo()
{
    int i, j;
    for(i = arinfocnt;i < infocnt;i++){
        for(j = 0;j < 3;j++){
            info[i][j].w = faceVector[j][infoToFaceIndex[i]];
            info[i][j + 3].w = -1 * faceVector[j][infoToFaceIndex[i]];
        }
    }
}
void L0ByFaceNeighborEdge::delInfo()
{
    int i;
    if(info){
        for(i =  0;i < infocnt;i++){
            if(info[i]){
                delete []info[i];
                info[i] = NULL;
            }
        }
        delete []info;
        info = NULL;
    }
}
