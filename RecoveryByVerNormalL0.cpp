#include "SubSolving.h"
#include "math.h"
#include "RecoveryByVerNormalL0.h"
#define pi acos(-1.0)

RecoveryByVerNormalL0::RecoveryByVerNormalL0(GLMmodel *pmeshmodel,IndexList **pverticesvindices,IndexList **pverticestindices,double **pvervector)
{
    meshmodel = pmeshmodel;
    verticestindices = pverticestindices;
    verticesvindices = pverticesvindices;
    vervector = pvervector;
}

RecoveryByVerNormalL0::~RecoveryByVerNormalL0()
{
    int i;
    if(relation)
    {
        for(i = 0;i < 3 * (int)meshmodel->numvertices;i++)
        {
            if(relation[i])
            {
                removeList(relation[i]);
            }
            relation[i] = NULL;
        }
        delete [] relation;
        relation = NULL;
    }

    if(info)
    {
        for(i = 0;i < (2 * edgeCnt) + (int)meshmodel->numvertices * 3;i++)
        {
            if(info[i])
                delete [] info[i];
            info[i] = NULL;
        }
        delete [] info;
        info = NULL;
    }
    if(p)
        delete []p;
    if(v)
        delete []v;

    relation = NULL;
    info = NULL;
    p = NULL;
    v = NULL;
}

void RecoveryByVerNormalL0::slove()
{
    int i, j;
    int maxtimes = 10;
    double beta = 0.7;
    double arpha = 1;
    double lambda = 0.0003;


    initRelation();
    initInfo();

    p = new double[3 * (int)meshmodel->numvertices + 2 * edgeCnt + arinfocnt];
    v = new double[3 * (int)meshmodel->numvertices];
    for(i = 0;i < 3 * (int)meshmodel->numvertices;i++)
    {
        p[i] = meshmodel->vertices[i + 3];
    }

    for(i = 0;i < (int)meshmodel->numvertices;i++){
        j = 0;
        p[i * 3 + j + (int)meshmodel->numvertices * 3] = p[3 * i + 1] * (-1) * vervector[j][2] + p[3 * i + 2] * vervector[j][1];

        j++;
        p[i * 3 + j + (int)meshmodel->numvertices * 3] = p[3 * i + 0] * vervector[j][2] + p[3 * i + 2] * (-1) * vervector[j][0];

        j++;
        p[i * 3 + j + (int)meshmodel->numvertices * 3] = p[3 * i + 0] * (-1) * vervector[j][1] + p[3 * i + 1] * vervector[j][0];
    }

    for(i = 0;i < 3 * (int)meshmodel->numvertices;i++)
    {
        v[i] = meshmodel->vertices[i + 3];
    }

    SubSolving s_l0(relation, info, 2 *  edgeCnt + arinfocnt, arinfocnt ,  3 * (int)meshmodel->numvertices);
    s_l0.init();

    int cc = 1;
    while(cc <= maxtimes)
    {
        for(i = arinfocnt;i < 2 * edgeCnt + arinfocnt;i++)
        {
            double sum = 0;
            for(j = 0;j < 6;j++)
            {
                sum += info[i][j].w * v[info[i][j].data];
            }
            if(sum * sum <= lambda / beta)
            {
                p[i + (int)meshmodel->numvertices * 3] = 0.0;
            }
            else
            {
                p[i + (int)meshmodel->numvertices * 3] = sum * beta;
            }
        }

        s_l0.getParameter(p, v,info, beta, arpha);
        s_l0.update();
        s_l0.getParameter(p, v,info, beta, arpha);
        s_l0.slove();
        printf("\trecoveryByL0\t%d\ttime\tfinished beta = %lf\n",cc,beta);
        beta = sqrt(2) * beta;
        cc++;
    }
    printf("\n\n");

	int yy[3] = {0};
    double xsum = -1;
    int idx = -1;

    for(i = 0;i < (int)meshmodel->numvertices;i++)
    {
        double len = 0;
        for(j = 0;j < 3;j++){
            len += (meshmodel->vertices[3 * i + j + 3] - v[3 * i + j]) * (meshmodel->vertices[3 * i + j + 3] - v[3 * i + j]);
        }
        len = sqrt(len);
        if(xsum < len){
            xsum = len;
            idx = i;
        }
    }
    double ysum  = 1.0;
    double xy = 0.0;
    xsum = xsum > 1e-3 ? xsum :1e-3;

    for(i = 0;i < 3;i++){
        xy += meshmodel->vertices[3 * idx + j + 3] * v[3 * idx + j];
    }
    double xxx = xy / ysum / xsum / 2.0;
    xxx = xxx > 1.0 ? 1.0 : xxx;
    xxx = xxx < -1.0 ? -1.0 : xxx;

    double angle = acos(xxx);

    if(angle > pi / 2){
        angle = pi - angle;
    }

    printf("max_diff : xsum = %lf angle = %lf\n",xsum, angle);

    for(i = 0;i < (int)meshmodel->numvertices;i++){
        for(j = 0;j < 3;j++){
            meshmodel->vertices[3 * i + j + 3] = v[3 * i + j];
        }
    }
}

void RecoveryByVerNormalL0::getParameter(List **pVerRelation,int pEdgeCnt, Edge *pEdge)
{
    verRelation = pVerRelation;
    edgeCnt = pEdgeCnt;
    edge = pEdge;

}

void RecoveryByVerNormalL0::initRelation()
{
    int i, j, k;
    relation = new List*[(int)meshmodel->numvertices * 3];
    for(i = 0;i < (int)meshmodel->numvertices * 3;i++)
    {
        relation[i] = NULL;
    }

    List *verrelationtail = NULL;
    List *relationtail[3];

    for(i = 0;i < (int)meshmodel->numvertices;i++)
    {
        for(j = 0;j < 3;j++)
        {
            verrelationtail = verRelation[i];
            while(verrelationtail)
            {
                for(k = 0;k < 3;k++)
                {
                    List *s = new List();
                    s->next = NULL;
                    s->data = 3 * verrelationtail->data + k;

                    if(NULL == relation[3 * i + j])
                    {
                        relation[3 * i + j] = s;
                        relationtail[j] = relation[3 * i + j];
                    }
                    else
                    {
                        relationtail[j]->next = s;
                        relationtail[j] = relationtail[j]->next;
                    }
                }
                verrelationtail = verrelationtail->next;
            }
        }
    }

    verrelationtail = NULL;
    for(i = 0;i < 3;i++)
        relationtail[i] = NULL;
}

void RecoveryByVerNormalL0::initInfo()
{
    int i, j;
    info = new Info *[(2 * edgeCnt) + (int)meshmodel->numvertices * 3];
    arinfocnt = (int)meshmodel->numvertices * 3;
    Edge *tail = edge;

    for(i = 0;i < (int)meshmodel->numvertices;i++){
        j = 0;
        info[i * 3 + j] = new Info[2];
        info[i * 3 + j][0].data = 3 * i + 1;
        info[i * 3 + j][0].w = -1 * vervector[j][2];
        info[i * 3 + j][0].cnt = 2;
        info[3 * i + j][1].data = 3 * i + 2;
        info[i * 3 + j][1].w = vervector[j][1];
        info[i * 3 + j][1].cnt = 2;

        j++;
        info[i * 3 + j] = new Info[2];
        info[i * 3 + j][0].data = 3 * i + 0;
        info[i * 3 + j][0].w = vervector[j][2];
        info[i * 3 + j][0].cnt = 2;
        info[i * 3 + j][1].data = 3 * i + 2;
        info[i * 3 + j][1].w = -1 * vervector[j][0];
        info[i * 3 + j][1].cnt = 2;

        j++;
        info[i * 3 + j] = new Info[2];
        info[i * 3 + j][0].data = 3 * i + 0;
        info[i * 3 + j][0].w = -1 * vervector[j][1];
        info[i * 3 + j][0].cnt = 2;
        info[i * 3 + j][1].data = 3 * i + 1;
        info[i * 3 + j][1].w = vervector[j][0];
        info[i * 3 + j][1].cnt = 2;
    }

    i = (int)meshmodel->numvertices * 3;
    while(tail)
    {
        info[i] = new Info[6];
		info[i + edgeCnt] = new Info[6];

        for(j = 0;j < 3;j++)
        {
            info[i][j].data = 3 * tail->p[0] + j;
            info[i][j].w = vervector[j][tail->p[0]];
            info[i][j].cnt = 6;

            info[i + edgeCnt][j].data = 3 * tail->p[0] + j;
            info[i + edgeCnt][j].w = vervector[j][tail->p[2]];
            info[i + edgeCnt][j].cnt = 6;
        }
        for(j = 0;j < 3;j++)
        {
            info[i][j + 3].data = 3 * tail->p[2] + j;
            info[i][j + 3].w = -1 * vervector[j][tail->p[0]];
            info[i][j + 3].cnt = 6;

            info[i + edgeCnt][j + 3].data = 3 * tail->p[2] + j;
            info[i + edgeCnt][j + 3].w = -1 * vervector[j][tail->p[2]];
            info[i + edgeCnt][j + 3].cnt = 6;
        }
        tail = tail->next;
        i++;
    }
    tail = NULL;
}
