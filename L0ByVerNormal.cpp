#include "L0ByVerNormal.h"
#include "RecoveryByVerNormalL0.h"

L0ByVerNormal::L0ByVerNormal(GLMmodel *originmeshmodel, IndexList **originverticesvindices, IndexList **originverticestindices) \
        :L0(originmeshmodel, originverticesvindices, originverticestindices)
{
    info = NULL;
    p = NULL;
    v = NULL;
}

L0ByVerNormal::~L0ByVerNormal()
{
    int i;
    if(p){
        for(i = 0;i < 3;i++){
            delete []p[i];
            p[i] = NULL;
        }
        delete []p;
        p = NULL;
    }

    if(v){
        for(i = 0;i < 3;i++){
            delete []v[i];
            v[i] = NULL;
        }
        delete []v;
        v = NULL;
    }
    delInfo();
    delVerNeighborVer();
    delVerSpreadNeighborVer();
}

GLMmodel* L0ByVerNormal::doL0(double parpha, double pbeta, double plambda, int pmaxtimes)
{

    int i,j,k;
    arpha = parpha;
    beta = pbeta;
    lambda = plambda;
    maxtimes = pmaxtimes;
    double tsum[3],sum;
    int cc = 1;

    //点的一环领域
    initVerNeighborVer();
    //点的二环领域
    initVerSpreadNeighborVer();
    //点法向
    initVerVector();
    updateVerVectorByArea();

    initInfo();
    updateInfo();
    getPV();

    SubSolving s_l0(verSpreadNeighborVer, info, infocnt, arinfocnt, (int)meshmodel->numvertices);
    s_l0.init();

    int nn, mm;
    while(cc <= pmaxtimes){
        mm = nn = 0;
        updateInfo();

        for(i = 0;i < (int)meshmodel->numvertices;i++){
            sum = 0.0;
            for(j = 0;j < 3;j++){
                tsum[j] = 0.0;
            }

            for(j = 0;j < info[i][j].cnt;j++){
                for(k = 0;k < 3;k++){
                    tsum[k] += verVector[k][info[i][j].data] * info[i][j].w;
                }
            }
            for(j = 0;j < 3;j++){
                sum += tsum[k] * tsum[k];
            }
            if(sum <= lambda / beta){
                nn++;
                for(k = 0;k < 3;k++)
                    p[k][i + (int)meshmodel->numvertices] = 0.0;
            }else{
                mm++;
                for(k = 0;k < 3;k++)
                    p[k][i + (int)meshmodel->numvertices] = beta * tsum[k];
            }
        }

        s_l0.getParameter(p[0], v[0], info, beta, arpha);
        s_l0.update();
        /**********************************************************/
        for(i = 0;i < 3;i++){
            //            传入p,v

            s_l0.getParameter(p[i], v[i], info, beta, arpha);
            s_l0.slove();
            for(j = 0;j < (int)meshmodel->numvertices;j++){
                verVector[i][j] = v[i][j];
            }
        }
        for(i = 0; i < (int)meshmodel->numvertices;i++){
            sum = 0.0;
            for(j = 0;j < 3;j++){
                sum += verVector[j][i] * verVector[j][i];
            }
            sum = sqrt(sum);
            sum = sum > 1e-12 ? sum : 1e-12;
            for(j = 0;j < 3;j++){
                verVector[j][i] = verVector[j][i] / sum;
            }
        }

		recoveryVerticesByVerNormal();
		printf("%d\ttime finished\n",cc);
        cc++;
        beta *= sqrt(2);
        arpha *= 2;
    }

    printf("finished\n");
	
    return meshmodel;
}

void L0ByVerNormal::getPV()
{
    int i, j;
    p = new double*[3];
    v = new double*[3];
    if(NULL == p || NULL == v){
        printf("error in L0ByVerNormal::getPV\n");
        return ;
    }
    for(i = 0;i < 3;i++){
        p[i] = new double[(int)meshmodel->numvertices + (int)meshmodel->numvertices];
        v[i] = new double[(int)meshmodel->numvertices];

        if(NULL == p[i] || NULL == v[i]){
            printf("error in L0ByVerNormal::getPV\n");
            return ;
        }
    }
    for(i = 0;i < 3;i++){
        for(j = 0;j < (int)meshmodel->numvertices;j++){
            p[i][j] = verVector[i][j];
            v[i][j] = 0.0;
        }
    }
    for(i = 0;i < 3;i++){
        for(j = 0;j < infocnt;j++){
            p[i][j + (int)meshmodel->numvertices] = 0.0;
        }
    }
}

void L0ByVerNormal::initInfo()
{
    List *verneighborvertail = NULL;
    int i, j;
    vector<int>v;

    infocnt = (int)meshmodel->numvertices;
    arinfocnt = 0;

    info = new Info*[infocnt];
    if(NULL == info){
        printf("error in L0ByVerNormal::initInfo\n");
        return ;
    }
    if(NULL == verNeighborVer){
        printf("error in L0ByVerNormal::initInfo  verNeighborVer\n");
        return ;
    }

    for(i = 0;i < (int)meshmodel->numvertices;i++){
        v.clear();
        verneighborvertail = verNeighborVer[i];
        while(verneighborvertail){
            v.push_back(verneighborvertail->data);
            verneighborvertail = verneighborvertail->next;
        }
        info[i] = new Info[(int)v.size()];
        if(NULL == info[i]){
            printf("error in L0ByVerNormal::info[%d\n]",i);
            return ;
        }
        for(j = 0;j < (int)v.size();j++){
            info[i][j].w = 0.0;
            info[i][j].data = v[j];
            info[i][j].cnt = (int)v.size();
        }
    }
/*    printf("info================\n");*/
    //for(i = 0; i< (int)meshmodel->numvertices;i++){
        //printf("%d\n",i);
        //for(j = 0;j < (int)info[i][j].cnt;j++){
            //printf("\t%lf %d %d\n",info[i][j].w,info[i][j].data,info[i][j].cnt);
        //}
    //}
    /*printf("----------------------\n");*/
    verneighborvertail = NULL;
}

void L0ByVerNormal::updateInfo()
{
    int i, j, k;
    vector<double>distv;
    vector<double>anglev;
    double sum = 0.0,tm = 0.0;

    for(i = 0;i < (int)meshmodel->numvertices;i++){
        distv.clear();
        anglev.clear();
        sum = 0.0;

        for(j = 0;j < info[i][j].cnt;j++){
            if(i != info[i][j].data){
                tm = getPointDistance(i, info[i][j].data);
                distv.push_back(tm);
                sum += tm;
            }
        }
        sum = sum > 1e-12 ? sum : 1e-12;

        k = 0;
        for(j = 0;j < info[i][j].cnt;j++){
            if(i != info[i][j].data){
 //               info[i][j].w = -1 * distv[k] / sum;
                info[i][j].w = -1 / (int)distv.size();
                k++;
            }else{
                info[i][j].w = 1;
            }
        }
    }
}

void L0ByVerNormal::delInfo()
{
    if(NULL != info){
        for(int i = 0; i < (int)meshmodel->numvertices;i++){
            delete []info[i];
            info[i] = NULL;
        }
        delete []info;
        info = NULL;
    }
}

void L0ByVerNormal::recoveryVerticesByVerNormal()
{
    printf("vernum = %d\ntrinum = %d\n",(int)meshmodel->numvertices,(int)meshmodel->numtriangles);
    RecoveryByVerNormalL0 myrecoverybyl0(meshmodel, verticesvindices, verticestindices, verVector);

    getEdge();
    initVerNeighborVer();

    myrecoverybyl0.getParameter(verNeighborVer, edgeNum, edge);
    myrecoverybyl0.slove();

    delVerNeighborVer();
    delEdge();
}
