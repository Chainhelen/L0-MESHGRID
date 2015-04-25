#include "SubSolving.h"

SubSolving::SubSolving(List **prelation,Info **pinfo,int pinfocnt,int parinfocnt, int prow)
{
    relation = prelation;
    info = pinfo;
    infocnt = pinfocnt;
    arinfocnt = parinfocnt;
    row = prow;
    AT = NULL;
    C = NULL;
}

SubSolving::~SubSolving()
{
    int i;
    relation = NULL;
    info = NULL;
    p = NULL;
    v = NULL;
    if(pb)
    {
        delete []pb;
    }
    pb = NULL;

    for(i = 0;i < row;i++)
    {
        removeList(AT[i]);
        removeList(C[i]);
    }
    if(AT)
        delete []AT;
    if(C)
        delete []C;
    AT = NULL;
    C = NULL;
}
void SubSolving::removeSparseMatrix(SparseMatrix *Head)
{
    if(NULL == Head)
        return;
    SparseMatrix *pNode = Head;
    SparseMatrix *pDel = NULL;

    while(pNode->next)
    {
        pDel = pNode->next;
        pNode->next = pDel->next;
        delete(pDel);
    }
    delete pNode;
    pNode = NULL;
    pDel = NULL;
}

void SubSolving::slove()
{
    mkl.MKLInitialization(C, row);
    AT[0]->ConstructNormalBvector(AT, p, pb, row);
    mkl.SolveSparseSystemMKL(pb, row, v);
}

void SubSolving::init()
{
    pb = new double[row];
    getATMem();
    getCMem();
}

void SubSolving::update()
{
    updateAT();
    updateC();
}

void SubSolving::getParameter(double *pp,double *pv,Info ** pinfo,double pbeta,double parpha)
{
    beta = pbeta;
    arpha = parpha;
    info = pinfo;
    p = pp;
    v = pv;
}

void SubSolving::getATMem()
{
    int i, j;
    SparseMatrix *sparsetail = NULL;

    AT = new SparseMatrix *[row];
    for(i = 0;i < row;i++)
    {
        SparseMatrix *s = new SparseMatrix();
        s->column = i;
        s->data = 1.0;
        s->next = NULL;
        AT[i] = s;
    }

    for(i = 0;i < infocnt;i++)
    {
        for(j = 0;j < info[i][0].cnt;j++)
        {
            sparsetail = AT[info[i][j].data];
            while(sparsetail->next)
            {
                sparsetail = sparsetail->next;
            }
            SparseMatrix *s = new SparseMatrix();
            s->column = i + row;
            s->next = NULL;
            s->data = 1.0;

            sparsetail->next = s;

            s = NULL;
        }
    }

    sparsetail = NULL;
}

void SubSolving::getCMem()
{
    int i;
    SparseMatrix *sparsetail = NULL;
    List *relationtail = NULL;

    C = new SparseMatrix *[row];
    for(i = 0;i < row;i++)
    {
        C[i] = NULL;
    }

    for(i = 0;i < row;i++)
    {
        relationtail = relation[i];
        while(relationtail)
        {
            SparseMatrix *s = new SparseMatrix();
            s->data = 1.0;
            s->column = relationtail->data;
            s->next = NULL;

            if(C[i] == NULL)
            {
                C[i] = s;
                sparsetail = C[i];
            }
            else
            {
                sparsetail->next = s;
                sparsetail = sparsetail->next;
            }

            s = NULL;
            relationtail = relationtail->next;
        }
    }
    sparsetail = NULL;
}

void SubSolving::updateAT()
{
    int i, j;
    SparseMatrix *sparsetail = NULL;

    for(i = 0;i < row;i++)
    {
        sparsetail = AT[i]->next;
        while(sparsetail)
        {
            for(j = 0;j < info[sparsetail->column - row][0].cnt;j++)
            {
                if(sparsetail->column - row - arinfocnt < 0)
                {
                    if(info[sparsetail->column - row][j].data == i)
                    {
                        sparsetail->data = info[sparsetail->column - row][j].w * arpha;
                        break;
                    }
                }
                if(sparsetail->column - row - arinfocnt >= 0)
                {
                    if(info[sparsetail->column - row][j].data == i)
                    {
                        sparsetail->data = info[sparsetail->column - row][j].w * beta;
                        break;
                    }
                }
            }
            sparsetail = sparsetail->next;
        }
    }
    sparsetail = NULL;
}

void SubSolving::updateC()
{
    int i;
    List *relationtail = NULL;
    SparseMatrix *ctail = NULL;
    for(i = 0;i < row;i++)
    {
        ctail = C[i];
        relationtail = relation[i];
        while(relationtail)
        {
            ctail->data = Sparsesparse(AT[i], AT[relationtail->data]);

            ctail = ctail->next;
            relationtail = relationtail->next;
        }
    }
    relationtail = NULL;
    ctail = NULL;
}

double SubSolving::Sparsesparse(SparseMatrix *A,SparseMatrix *B)
{
	SparseMatrix *point = A;
	SparseMatrix *point2 = B;

	double sum = 0.0;

	while(point)
	{
		while(point2)
		{
			if(point->column == point2->column)
			{
				sum += point->data * point2->data;
			}
			else if(point->column < point2->column)
			{
				break;
			}
			point2 = point2->next;
		}
		point = point->next;
	}
	point = NULL;
	point2 = NULL;

	return sum;
}
