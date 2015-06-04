#include "SubSolving.h"
#include "math.h"
#include "RecoveryByVerNormalL0.h"

#define ADD_ARPHA_FIRST
//#define ADD_ARPHA_SECOND

RecoveryByVerNormalL0::RecoveryByVerNormalL0(GLMmodel *pmeshmodel,List **pverticesortvindices,double **pvervector)
{
    meshmodel = pmeshmodel;
    verticesortvindices = pverticesortvindices;
    vervector = pvervector;
    infocnt = 0;
    arinfocnt = 0;
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
        for(i = 0;i < edgeCnt;i++)
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
    double beta = 0.07;
    double arpha = 2;
    double lambda = 0.003;

    initRelation();
    initInfo();

    p = new double[3 * (int)meshmodel->numvertices + infocnt];
    v = new double[3 * (int)meshmodel->numvertices];
    for(i = 0;i < 3 * (int)meshmodel->numvertices;i++)
    {
        p[i] = meshmodel->vertices[i + 3];
    }
    for(i = 0; i < arinfocnt;i++){
        p[i + 3 * (int)meshmodel->numvertices] = 0.0;
    }
    for(i = 0;i < 3 * (int)meshmodel->numvertices;i++)
    {
        v[i] = meshmodel->vertices[i + 3];
    }

    SubSolving s_l0(relation, info, infocnt, arinfocnt ,  3 * (int)meshmodel->numvertices);
    s_l0.init();

    int cc = 1;
    while(cc <= maxtimes)
    {
        for(i = arinfocnt;i < infocnt;i++)
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
        printf("\trecoveryByL0\t%d\ttime\tfinished\n",cc);
//		arpha /= 2;
        beta = sqrt(2) * beta;
        cc++;
    }
    printf("\n\n");
    for(i = 0;i < 3 * (int)meshmodel->numvertices;i++)
    {
        meshmodel->vertices[i + 3] = v[i];
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
    arinfocnt = 0;
    Edge *tail = NULL;

#ifdef ADD_ARPHA_FIRST
    arinfocnt = edgeCnt;
#endif

#ifdef ADD_ARPHA_SECOND
    arinfocnt = 3 * meshmodel->numvertices;
#endif

    infocnt = 2 * edgeCnt + arinfocnt;
    info = new Info *[infocnt];


#ifdef ADD_ARPHA_FIRST
    tail = edge;
    i = 0;
    while(tail){
        info[i] = new Info[12];
        for(j = 0;j < 4;j++){
            for(int k = 0;k < 3;k++){
                info[i][3 * j + k].cnt = 12;
                info[i][3 * j + k].w = (j % 2 == 0 ? -1 : 1);
                info[i][3 * j + k].data = 3 * tail->p[j] + k;
            }
        }
        i++;
        tail = tail->next;
    }
#endif

#ifdef ADD_ARPHA_SECOND
    for(i = 0;i < meshmodel->numvertices;i++){
        List *vertail = verticesortvindices[i];
        int num = 0;

        while(vertail){
            num++;
            vertail = vertail->next;
        }

        for(j = 0;j < 3;j++){
            info[3 * i + j] = new Info[num];
        }
        vertail = verticesortvindices[i];

        int k = 0;
        while(vertail){
            if(vertail->data == i){
                for(j = 0;j < 3;j++){
                    info[3 * i + j][k].data = 3 *vertail->data + j;
                    info[3 * i + j][k].w = -1.0;
                    info[3 * i + j][k].cnt = num;
                }
            }else{
                for(j = 0;j < 3;j++){
                    info[3 * i + j][k].data = 3 * vertail->data + j;
                    info[3 * i + j][k].w = 1.0 / (num - 1);
                    info[3 * i + j][k].cnt = num;
                }
            }
            k++;
            vertail = vertail->next;
        }
    }
#endif

    tail = edge;
    i = arinfocnt;
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
