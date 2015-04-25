#include "L0BYEDGE.H"

L0ByEdge::L0ByEdge(GLMmodel *originmeshmodel, IndexList **originverticesvindices, IndexList **originverticestindices) \
        :L0(originmeshmodel, originverticesvindices, originverticestindices)
{
    info = NULL;
    p = NULL;
    v = NULL;
}

L0ByEdge::~L0ByEdge()
{
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
}

GLMmodel* L0ByEdge::doL0(double parpha, double pbeta, double plambda, int pmaxtimes)
{
    arpha = parpha;
    beta = pbeta;
    lambda = plambda;
    maxtimes = pmaxtimes;
    int cc = 1, i, j, k;
    double tsum[3], sum;
    //得到边,求info
    getEdge();
    initInfo();
    delEdge();

    getPV();
    //得到排序后二环领域，包括自身点
    initVerSpreadNeighborVer();

    //稀疏求解初始化
    SubSolving s_l0(verSpreadNeighborVer, info, infocnt, arinfocnt, (int)meshmodel->numvertices);
    s_l0.init();

    while(cc <= maxtimes){
        updateInfo();

        for(i = arinfocnt;i < infocnt;i++){
            sum = 0;
            for(j = 0;j < 3;j++)
                tsum[j] = 0;

            for(j = 0;j < 4;j++)
            {
                for(k = 0;k < 3;k++){
                    tsum[k] += info[i][j].w * model[info[i][j].data * 3 + k];
                }
            }
            for(j = 0;j < 3;j++){
                sum += tsum[j] * tsum[j];
            }
            if(sum <= lambda / beta){
                for(j = 0;j < 3;j++)
                    p[j][i + (int)meshmodel->numvertices] = 0.0;
            }
            else{
                for(j = 0;j < 3;j++)
                    p[j][i + (int)meshmodel->numvertices] = beta * tsum[j];
            }
        }

        //传入beta,arpha
        s_l0.getParameter(p[0], v[0], info, beta, arpha);
        s_l0.update();

        for(i = 0;i < 3;i++){
            //传入p,v
            s_l0.getParameter(p[i], v[i], info, beta, arpha);
            s_l0.slove();
            for(j = 0;j < (int)meshmodel->numvertices;j++){
                model[3 * j + i] = v[i][j];
            }
        }

        printf("%d\t time finished\n",cc);
        cc++;
        beta = beta * sqrt(2);
        arpha = arpha / 2;
    }
    //删除info,和关系链表
    delInfo();
    delVerSpreadNeighborVer();

    for(i = 0;i < (int)meshmodel->numvertices;i++){
        for(j = 0;j < 3;j++){
            meshmodel->vertices[3 * (i + 1) + j]  = model[3 * i + j];
        }
    }
    return meshmodel;
}

void L0ByEdge::getPV()
{
    int i, j;
    p = new double*[3];
    v = new double*[3];
    if(NULL == p || NULL == v){
        printf("Error in L0ByEdge::doL0----- p , v\n");
        return ;
    }

    for(i = 0;i < 3;i++){

        p[i] = new double[(int)meshmodel->numvertices + infocnt];
        v[i] = new double[(int)meshmodel->numvertices];

        if(NULL == p[i] || NULL == v[i]){
            printf("Error in L0ByEdge::doL0----- p[%d] , v[%d]\n",i,i);
            return ;
        }
    }
    for(i = 0; i < 3;i++){
        for(j = 0;j < (int)meshmodel->numvertices;j++){
            p[i][j] = model[3 * j + i];
            v[i][j] = model[3 * j + i];
        }
    }
    for(i = (int)meshmodel->numvertices;i < (int)meshmodel->numvertices + arinfocnt;i++){
        for(j = 0;j < 3;j++){
            p[j][i] = 0;
        }
    }
}

void L0ByEdge::initInfo()
{
    int i, j;
    arinfocnt = edgeNum;
    infocnt = 2 * edgeNum;

    if(NULL != info){
        delInfo();
    }

    info = new Info*[infocnt];

    Edge *edgetail  = edge;

    if(NULL == edge || NULL == info){
        printf("error in L0ByEdge::getInfo\n");
        return ;
    }

    for(i = 0;i < edgeNum;i++){
        info[i] = new Info[4];

        if(NULL == info[i]){
            printf("error in L0ByEdge::getInfo----->info[i]\n");
            return ;
        }

        info[i + arinfocnt] = new Info[4];

        if(NULL == info[i]){
            printf("error in L0ByEdge::getInfo----->info[i + arinfocnt]\n");
            return ;
        }

        for(j = 0;j < 4;j++){
            info[i][j].data = edgetail->p[j];
            info[i][j].cnt = 4;
            info[i][j].w = 0.0;

            info[i + arinfocnt][j].data = edgetail->p[j];
            info[i + arinfocnt][j].cnt = 4;
            info[i + arinfocnt][j].w = 0.0;
        }
        //printf("%d \n",edgetail->p[0]);
        edgetail = edgetail->next;
    }
    edgetail = NULL;
}

void L0ByEdge::updateInfo()
{
    int i, j;
    for(i = 0;i < arinfocnt;i++)
    {
        for(j = 0;j < 4;j++)
        {
            if(j % 2 ==	0)
                info[i][j].w = 1 * arpha;
            else
                info[i][j].w = -1 * arpha;
        }
    }

    for(i = arinfocnt;i < infocnt;i++)
    {
        int p1 = info[i][0].data;
        int p2 = info[i][1].data;
        int p3 = info[i][2].data;
        int p4 = info[i][3].data;

        double p1p2 = getPointDistance(p1, p2);
        double p2p3 = getPointDistance(p2, p3);
        double p1p3 = getPointDistance(p1, p3);
        double p1p4 = getPointDistance(p1, p4);
        double p3p4 = getPointDistance(p3, p4);

        double angle_312 = getPointAngle(p2p3, p1p2, p1p3);
        double angle_231 = getPointAngle(p1p2, p2p3, p1p3);
        double angle_413 = getPointAngle(p3p4, p1p3, p1p4);
        double angle_134 = getPointAngle(p1p4, p1p3, p3p4);

        double area123 = 0.5 * p2p3 * p1p3 * sin(angle_231);
        double area134 = 0.5 * p3p4 * p1p3 * sin(angle_134);

        double fenmu;
        fenmu = (p1p3) * (p1p3) * (area123 + area134);
        if(fenmu < 1e-14)
            fenmu = 1e-14;
        info[i][0].w = (area123 * getPointFactor(p4,p3,p3,p1) + area134 * getPointFactor(p1,p3,p3,p2)) / (fenmu);

        fenmu = area123 + area134;
        if(fenmu < 1e-14)
            fenmu = 1e-14;
        info[i][1].w = area134 / (fenmu);

        fenmu = (p1p3) * (p1p3) * (area123 + area134);
        if(fenmu < 1e-14)
            fenmu = 1e-14;
        info[i][2].w = (area123 * getPointFactor(p3,p1,p1,p4) + area134 * getPointFactor(p2,p1,p1,p3)) / (fenmu);

        fenmu = area123 + area134;
        if(fenmu < 1e-14)
            fenmu = 1e-14;
        info[i][3].w = area123 / (fenmu);
    }
}

void L0ByEdge::delInfo()
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
