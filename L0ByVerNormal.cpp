#include "L0ByVerNormal.h"
#include "RecoveryByVerNormalL0.h"

L0ByVerNormal::L0ByVerNormal(GLMmodel *originmeshmodel, IndexList **originverticesvindices, IndexList **originverticestindices) \
        :L0(originmeshmodel, originverticesvindices, originverticestindices)
{
    info = NULL;
    p = NULL;
    v = NULL;
    guidedmesh = NULL;
    guidedmeshVerNormal = NULL;
    guidedmeshFaceNormal = NULL;
}

void L0ByVerNormal::getGuidedmesh(GLMmodel *pGuidedmesh){
    guidedmesh = pGuidedmesh;
}

void L0ByVerNormal::initGuidedmeshFaceNormal(){
    int tindex, i, j;
    double vec[3];

    guidedmeshFaceNormal = new double *[3];
    for(i = 0;i < 3;i++){
        guidedmeshFaceNormal[i] = new double [(int)meshmodel->numtriangles];
    }
    printf("%lf %lf %lf\n",guidedmesh->vertices[3],guidedmesh->vertices[4],guidedmesh->vertices[5]);

    for(tindex = 0;tindex < (int)meshmodel->numtriangles;tindex++)
    {
        double u[3], v[3];

        int a = meshmodel->triangles[tindex].vindices[0];
        int b = meshmodel->triangles[tindex].vindices[1];
        int c = meshmodel->triangles[tindex].vindices[2];

        u[0] = guidedmesh->vertices[3 * b + 0] - guidedmesh->vertices[3 * a + 0];
        u[1] = guidedmesh->vertices[3 * b + 1] - guidedmesh->vertices[3 * a + 1];
        u[2] = guidedmesh->vertices[3 * b + 2] - guidedmesh->vertices[3 * a + 2];

        v[0] = guidedmesh->vertices[3 * c + 0] - guidedmesh->vertices[3 * a + 0];
        v[1] = guidedmesh->vertices[3 * c + 1] - guidedmesh->vertices[3 * a + 1];
        v[2] = guidedmesh->vertices[3 * c + 2] - guidedmesh->vertices[3 * a + 2];


        vec[0] = u[1] * v[2] - u[2] * v[1];
        vec[1] = u[2] * v[0] - u[0] * v[2];
        vec[2] = u[0] * v[1] - u[1] * v[0];


        double sum = 0.0;
        for(j = 0;j < 3;j++){
            sum += vec[j] * vec[j];
        }
        sum = sqrt(sum);
        sum = sum > 1e-12 ? sum : 1e-12;
        for(j = 0;j < 3;j++){
            vec[j] /= sum;
        }
        for(j = 0;j < 3;j++){
            guidedmeshFaceNormal[j][tindex] = vec[j];
        }
    }
}

void L0ByVerNormal::initGuidedmeshVerNormal(){
    vector<int>v;
    double tsum[3];
    double sum;
    int i, j;


    guidedmeshVerNormal = new double *[3];
    for(i = 0;i < 3;i++){
        guidedmeshVerNormal[i] = new double [(int)meshmodel->numvertices];
    }

    for(i = 0;i < (int)meshmodel->numvertices;i++){
        v.clear();
        IndexList *tail = verticestindices[i + 1];

        for(j = 0;j < 3;j++){
            tsum[j] = 0.0;
            sum = 0.0;
        }

        while(tail){
            for(j = 0;j < 3;j++){
                tsum[j] += guidedmeshFaceNormal[j][tail->index];
            }
            tail = tail->next;
        }

        for(j = 0;j < 3;j++){
            sum += tsum[j] * tsum[j];
        }
        sum = sqrt(sum);
        sum = sum > 1e-12 ? sum : 1e-12;

        for(j = 0;j < 3;j++){
            guidedmeshVerNormal[j][i] = tsum[j] / sum;
        }
    }
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
    initGuidedmeshFaceNormal();
    initGuidedmeshVerNormal();

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
    nn = 0;
    mm = 0;
    for(i = 0;i < (int)meshmodel->numvertices;i++){
        double sum = 0;
        for(j = 0;j < 3;j++){
            sum += verVector[j][i] * verVector[j][i];
        }
        sum = sqrt(sum);
        if(fabs(sum - 1.0) > 1e-6){
            nn++;
            printf("find1 \n");
        }
        sum = 0;
        for(j = 0;j < 3;j++){
            sum += guidedmeshVerNormal[j][i] * guidedmeshVerNormal[j][i] ;
        }
        sum = sqrt(sum);
        if(fabs(sum - 1.0) > 1e-6){
            mm++;
            printf("find2 \n");
        }
    }
    printf("nn mm %d %d\n",nn,mm);


	recoveryVerticesByVerNormal();
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
    //printf("vernum = %d\ntrinum = %d\n",(int)meshmodel->numvertices,(int)meshmodel->numtriangles);
    RecoveryByVerNormalL0 myrecoverybyl0(meshmodel, verticesvindices, verticestindices, guidedmeshVerNormal);

    getEdge();
    initVerNeighborVer();

    myrecoverybyl0.getParameter(verNeighborVer, edgeNum, edge);
    myrecoverybyl0.slove();

    delVerNeighborVer();
    delEdge();
}