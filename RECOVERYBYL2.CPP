#include "RECOVERYBYL2.H"

RecoveryByL2::RecoveryByL2(GLMmodel *originmeshmodel, IndexList **originverticesvindices, IndexList **originverticestindices) \
        :L0(originmeshmodel, originverticesvindices, originverticestindices)
{
    info = NULL;
    verThreeSeparateRelation = NULL;
    p = NULL;
    v = NULL;
    model = NULL;
    faceVector = NULL;
}

RecoveryByL2::~RecoveryByL2()
{

    if(NULL != p){
        delete []p;
        p = NULL;
    }
    if(NULL != v){
        delete []v;
        v = NULL;
    }
    model = NULL;
    faceVector = NULL;

    delInfo();
    delVerThreeSeparateRelation();
}

void RecoveryByL2::changeModel(double *pModel)
{
    model = pModel;
}

void RecoveryByL2::changeFaceVector(double **pFaceVector)
{
    faceVector = pFaceVector;
}

void RecoveryByL2::doL2(double pbeta)
{
    int i, j;
    beta = pbeta;
    p = new double[(int)meshmodel->numvertices * 3 + (int)meshmodel->numtriangles * 3];
    v = new double[(int)meshmodel->numvertices * 3];

    initInfo();
    updateInfo();

    initVerThreeSeparateRelation();

    for(i = 0;i < (int)meshmodel->numvertices;i++){
        for(j = 0;j < 3;j++){
            p[i * 3 + j] = model[3 * i + j];
        }
    }
    for(i = (int)meshmodel->numvertices * 3;i < (int)meshmodel->numvertices * 3 + (int)meshmodel->numtriangles * 3;i++){
        p[i] = 0;
    }
    SubSolving s_l0(verThreeSeparateRelation, info, infocnt, arinfocnt, (int)meshmodel->numvertices * 3);
    s_l0.init();

    //传入beta,arpha
    s_l0.getParameter(p, v, info, beta, arpha);
    s_l0.update();
    s_l0.slove();

    for(i = 0;i < (int)meshmodel->numvertices;i++){
        for(j = 0; j < 3;j++)
            model[3 * i + j] = v[3 * i + j];
    }

    delInfo();
    delVerThreeSeparateRelation();
}


void RecoveryByL2::initInfo()
{
    int i, j, k;
    arinfocnt = 0;
    infocnt = (int)meshmodel->numtriangles * 3;
    int aIndex, bIndex;

    info = new Info *[infocnt];
    if(NULL == info){
        printf("error in RecoveryByL2::initInfo\n");
        return ;
    }

    for(i = 0;i < (int)meshmodel->numtriangles;i++){
        for(j = 0;j < 3;j++){
            aIndex = meshmodel->triangles[i].vindices[j % 3] - 1;
            bIndex = meshmodel->triangles[i].vindices[(j + 1) % 3] - 1;

            info[3 * i + j] = new Info[6];
            if(NULL == info[3 * i + j]){
                printf("error in RecoveryByL2::initInfo info[%d]\n",3 * i + j);
            }
            for(k = 0;k < 3;k++){
                info[3 * i + j][k].data = 3 * aIndex + k;
                info[3 * i + j][k].cnt = 6;
                info[3 * i + j][k].w = 0.0;

                info[3 * i + j][k + 3].data = 3 * bIndex + k;
                info[3 * i + j][k + 3].cnt = 6;
                info[3 * i + j][k + 3].w = 0.0;
            }
        }
    }
}
void RecoveryByL2::updateInfo()
{
    int i, j, k;
    for(i = 0;i < (int)meshmodel->numtriangles;i++){
        for(j = 0;j < 3;j++){

            for(k = 0;k < 3;k++){
                info[3 * i + j][k].w = faceVector[k][i];
                info[3 * i + j][k + 3].w = -1 * faceVector[k][i];
            }
        }
    }
}

void RecoveryByL2::delInfo()
{
    if(info){
        for(int i = 0;i < infocnt;i++){
            delete []info[i];
            info[i] = NULL;
        }
        delete []info;
        info = NULL;
    }
}

void RecoveryByL2::initVerThreeSeparateRelation()
{

        List *verneighborvertail = NULL;
        int i, j, num;
        vector<int>v;
        verThreeSeparateRelation = new List*[(int)meshmodel->numvertices * 3];
        if(NULL == verThreeSeparateRelation){
            printf("error in RecoveryByL2::initVerThreeSeparateRelation ------->\n");
            return;
        }
        for(i = 0;i < (int)meshmodel->numvertices * 3;i++){
            verThreeSeparateRelation[i] = NULL;
        }

        initVerNeighborVer();

        for(i = 0;i < (int)meshmodel->numvertices;i++){
            verneighborvertail = verNeighborVer[i];
            v.clear();
            num = 0;

            while(verneighborvertail){
                for(j = 0;j < 3;j++){
                    v.push_back(3 * verneighborvertail->data + j);
                    num++;
                }
                verneighborvertail = verneighborvertail->next;
            }
            for(j = 0;j < 3;j++){
                insertList(verThreeSeparateRelation, 3 * i + j, v, num);
            }
        }

        delVerNeighborVer();
}

void RecoveryByL2::delVerThreeSeparateRelation()
{
    if(verThreeSeparateRelation){
        for(int i = 0;i < (int)meshmodel->numtriangles;i++){
            if(NULL != verThreeSeparateRelation[i]){
                removeList(verThreeSeparateRelation[i]);
            }
            verThreeSeparateRelation[i] = NULL;
        }
        delete []verThreeSeparateRelation;
        verThreeSeparateRelation = NULL;
    }
}
