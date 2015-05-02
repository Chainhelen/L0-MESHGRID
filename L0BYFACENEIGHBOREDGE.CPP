#include "L0BYFACENEIGHBOREDGE.H"

//#define FIRSTFACENEIGHBORFACE
#define SECONDFACENEIGHBORFACE
//#define ALLFACENEIGHBORFACE


L0ByFaceNeighborEdge::L0ByFaceNeighborEdge(GLMmodel *originmeshmodel, IndexList **originverticesvindices, IndexList **originverticestindices) \
        :L0(originmeshmodel, originverticesvindices, originverticestindices)
{
    info = NULL;
    faceNeighborFace = NULL;
    infoToFaceIndex = NULL;
    verNeighborThreeVer = NULL;
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

    delVector(faceVector);
    faceVector = NULL;

    delVerNeighborThreeVer();

    info = NULL;
    faceNeighborFace = NULL;
}

GLMmodel* L0ByFaceNeighborEdge::doL0(double parpha, double pbeta, double plambda, int pmaxtimes)
{
    int i, j;
    arpha = parpha;
    beta = pbeta;
    lambda = plambda;
    maxtimes = pmaxtimes;
    int cc = 1;
    double sum;
    List *tail = NULL;

    initVerNeighborVer();
    initVerNeighborThreeVer();

    chooseTypeOfNeighbor();
    initFaceVector();
    updateFaceVector();

    initInfo();
    updateInfo();

    getPV();

    SubSolving s_l0(verNeighborThreeVer, info, infocnt, arinfocnt, (int)meshmodel->numvertices * 3);
    s_l0.init();

    while(cc <= pmaxtimes){
        updateFaceVector();
        updateInfo();
        for(i = arinfocnt;i < infocnt;i++){
            sum = 0.0;
            for(j = 0;j < 6;j++){
                sum += model[info[i][j].data] * info[i][j].w;
            }
            if(sum * sum <= lambda / beta){
                p[i + (int)meshmodel->numvertices * 3] = 0.0;
            }else{
                p[i + (int)meshmodel->numvertices * 3] = sum * beta;
            }
        }
        s_l0.getParameter(p, v, info, beta, arpha);
        s_l0.update();
        s_l0.slove();

        for(i = 0;i < (int)meshmodel->numvertices;i++){
            for(j = 0;j < 3;j++){
                model[3 * i + j] = v[3 * i + j];
            }
        }
		printf("%d time\n",cc);
        cc++;
        beta = beta * sqrt(2);
        arpha = arpha / 2;
    }

    for(i = 0;i < (int)meshmodel->numvertices;i++){
        for(j = 0;j < 3;j++){
            meshmodel->vertices[3 * (i + 1) + j] = model[3 * i + j];
        }
    }
    return meshmodel;
}

void L0ByFaceNeighborEdge::initVerNeighborThreeVer()
{
    List *verneighborvertail = NULL;
    int i, j, num;
    vector<int>v;
    verNeighborThreeVer = new List*[(int)meshmodel->numvertices * 3];

    if(NULL == verNeighborThreeVer){
        printf("L0ByFaceNeighborEdge::initVerNeighborThreeVer\n");
        return;
    }

    for(i = 0;i < (int)meshmodel->numvertices * 3;i++){
        verNeighborThreeVer[i] = NULL;
    }

    for(i = 0;i < (int)meshmodel->numvertices;i++){
        verneighborvertail = verNeighborVer[i];
        v.clear();
        while(verneighborvertail){
            for(j = 0;j < 3;j++){
                v.push_back(3 * verneighborvertail->data + j);
            }
            verneighborvertail =verneighborvertail->next;
        }
        num  = (int)v.size();
        for(j = 0;j < 3;j++){
            insertList(verNeighborThreeVer, 3 * i + j, v, num);
        }
    }
    verneighborvertail = NULL;
}

void L0ByFaceNeighborEdge::delVerNeighborThreeVer()
{
    int i;

	if(verNeighborThreeVer){
		for(i = 0;i < (int)meshmodel->numvertices * 3;i++){
			if(verNeighborThreeVer[i])
				removeList(verNeighborThreeVer[i]);
			verNeighborThreeVer[i] = NULL;
		}
	}
    if(verNeighborThreeVer){
        delete []verNeighborThreeVer;
        verNeighborThreeVer = NULL;
    }
}

void L0ByFaceNeighborEdge::getPV()
{
    p = new double [(int)meshmodel->numvertices * 3 + infocnt];
    v = new double [(int)meshmodel->numvertices * 3];

    if(NULL == p || NULL == v){
        printf("error at L0ByFaceNeighborEdge::getPV\n");
        return ;
    }
    int i, j;

    for(i = 0;i < (int)meshmodel->numvertices;i++){
		for(j = 0;j < 3;j++)
		{
			p[3 * i + j] = model[3 * i + j];
			v[3 * i + j] = model[3 * i + j];
		}
    }
    for(i = 0;i < infocnt;i++){
        p[i + (int)meshmodel->numvertices * 3] = 0.0;
    }
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
                            if(b == -1){
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
                            if(b == -1){
                                b = meshmodel->triangles[fIndex].vindices[j];
                            }else{
                                c = meshmodel->triangles[fIndex].vindices[j];
                            }
                        }
                    }
                }
                a--,b--,c--;

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
    double len;
    for(i = arinfocnt;i < infocnt;i++){
        len = getPointDistance(info[i][0].data / 3, info[i][3].data / 3);
        if(len < 1e-12){
            printf("ll\n");
        }
        len = len > 1e-12 ? len : 1e-12;

        for(j = 0;j < 3;j++){
            info[i][j].w = faceVector[j][infoToFaceIndex[i]] / len;
            info[i][j + 3].w = -1 * faceVector[j][infoToFaceIndex[i]] / len;
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
