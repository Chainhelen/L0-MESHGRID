#include "L0BYFACEVECTOR.H"

//#define FIRSTFACENEIGHBORFACE
#define SECONDFACENEIGHBORFACE
//#define ALLFACENEIGHBORFACE

//#define TESTZHENG_LOCAL
//#define TEST2
#define  MAX_TMP 3

L0ByFaceVector::L0ByFaceVector(GLMmodel *originmeshmodel, IndexList **originverticesvindices, IndexList **originverticestindices) \
        :L0(originmeshmodel, originverticesvindices, originverticestindices)
{
    info = NULL;
    faceNeighborFace = NULL;
    faceSpreadNeighborFace = NULL;
    p = NULL;
    v = NULL;
}

L0ByFaceVector::~L0ByFaceVector()
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
    delFaceSpreadNeighborFace();

    faceNeighborFace = NULL;
    faceSpreadNeighborFace = NULL;
}

GLMmodel* L0ByFaceVector::doL0(double parpha, double pbeta, double plambda, int pmaxtimes)
{
    arpha = parpha;
    beta = pbeta;
    lambda = plambda;
    maxtimes = pmaxtimes;
    int cc = 1, i, j, k;
    double tsum[3], sum;


    chooseTypeOfNeighbor();

    initInfo();
    updateInfo();


    initFaceVector();
    updateFaceVector();

#ifdef TESTZHENG_LOCAL
    vector<double>distv;
    vector<double>normv;
    vector<double>areav;
    vector<double>wv;
    double w;

    double max_norm,max_dist;
    double distw, normw, areaw;
    double wsum;

    int ftemp = MAX_TMP;
    while(ftemp--){
        updateFaceVector();
        for(i = 0;i < (int)meshmodel->numtriangles;i++){
            distv.clear();
            normv.clear();
            areav.clear();
            max_dist = max_norm = -1;

            for(j = 0;j < 3;j++){
                tsum[j] = 0;
            }
            List *faceneighborfacetail = faceNeighborFace[i];
            while(faceneighborfacetail){
                if(i != faceneighborfacetail->data){
                    areaw = getAreaByIndex(faceneighborfacetail->data);
                    distw = getDistanceFaceCenterToFaceCenter(i, faceneighborfacetail->data);
                    normw =  getSignalByFaceIndex(i, faceneighborfacetail->data);

                    max_dist = max(max_dist, distw);
                    max_norm = max(max_norm, normw);

                    areav.push_back(areaw);
                    distv.push_back(distw);
                    normv.push_back(normw);
                }
                faceneighborfacetail =  faceneighborfacetail->next;
            }
            max_dist = max_dist > 1e-12 ? max_dist : 1e-12;
            max_norm = max_norm > 1e-12 ? max_norm : 1e-12;

            wsum = 0.0;
            wv.clear();
            for(j = 0;j < (int)areav.size();j++){
                w = areav[j] * exp(- distv[j] * distv[j] / (2 * max_dist * max_dist) - normv[j] * normv[j] / (2 * max_norm * max_norm));
                wsum += w;
                wv.push_back(w);
            }

            for(j = 0;j < 3;j++){
                tsum[j] = 0.0;
            }
            faceneighborfacetail = faceNeighborFace[i];
            k = 0;
            while(faceneighborfacetail){
                if(i != faceneighborfacetail->data){
                    for(j = 0;j < 3;j++){
                        tsum[j] += wv[k] / wsum * faceVector[j][faceneighborfacetail->data];
                    }
                    k++;
                }
                faceneighborfacetail = faceneighborfacetail->next;
            }
            for(j = 0;j < 3;j++){
                faceVector[j][i] = tsum[j];
            }
        }
        recoveryVerticesByFaceVector();
        printf("%d\n",ftemp);
    }
    for(i = 0;i < (int)meshmodel->numvertices;i++){
        for(j = 0;j < 3;j++){
            meshmodel->vertices[3 * (i + 1) + j]  = model[3 * i + j];
        }
    }

    return meshmodel;

#endif

    p = new double*[3];
    v = new double*[3];

    if(NULL == p || NULL == v){
        return NULL;
    }

    for(i = 0;i < 3;i++){

        p[i] = new double[(int)meshmodel->numtriangles + infocnt];
        v[i] = new double[(int)meshmodel->numtriangles];

        if(NULL == p[i] || NULL == v[i]){
            printf("Error in L0ByEdge::doL0----- p[%d] , v[%d]\n",i,i);
            return NULL;
        }
    }

    for(i = 0;i < 3;i++){
        for(j = 0;j < (int)meshmodel->numtriangles;j++){
            p[i][j] = faceVector[i][j];
            v[i][j] = faceVector[i][j];
        }
    }

    for(i = 0;i < 3;i++){
        for(j = (int)meshmodel->numtriangles;j < (int)meshmodel->numtriangles + arinfocnt;j++){
            p[i][j] = 0.0;
        }
    }

    SubSolving s_l0(faceSpreadNeighborFace, info, infocnt, arinfocnt, (int)meshmodel->numtriangles);
    s_l0.init();

    int nn, mm;
    int tempp = MAX_TMP;
    while(tempp--){
/*        for(i = 0;i < (int)meshmodel->numvertices;i++){*/
            //for(j = 0;j < 3;j++){
                //p[j][i] = faceVector[j][i];
            //}
        /*}*/

       arpha = parpha;
       beta = pbeta;
       lambda = plambda;
       cc = 1;
       updateFaceVector();
       updateInfo();

        while(cc <= maxtimes){
            nn = 0;
            mm = 0;
            for(i = arinfocnt;i < infocnt;i++){
                 sum = 0;
                 for(j = 0;j < 3;j++)
                     tsum[j] = 0;

                 for(j = 0;j < info[i][0].cnt;j++)
                 {
                     for(k = 0;k < 3;k++){
                         tsum[k] += info[i][j].w * faceVector[k][info[i][j].data];
                     }
                 }
                 for(j = 0;j < 3;j++){
                     sum += tsum[j] * tsum[j];
                 }
                 if(sum <= lambda / beta){
                     nn++;
                     for(j = 0;j < 3;j++)
                         p[j][i + (int)meshmodel->numtriangles] = 0.0;
                 }
                 else{
                     mm++;
                     for(j = 0;j < 3;j++)
                         p[j][i + (int)meshmodel->numtriangles] = beta * tsum[j];
                 }
            }
               //        传入beta,arpha

            s_l0.getParameter(p[0], v[0], info, beta, arpha);
            s_l0.update();
            /**********************************************************/
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
			recoveryVerticesByFaceVector();
            printf("%d %d\n",nn,mm);
            printf("%d\t time finished\n",cc);
            cc++;
            beta = beta * sqrt(2);
            arpha = arpha / 2;
        }
    }


        //   删除info,和关系链表
    delInfo();
    delFaceNeighborFace();

    for(i = 0;i < (int)meshmodel->numvertices;i++){
        for(j = 0;j < 3;j++){
            meshmodel->vertices[3 * (i + 1) + j]  = model[3 * i + j];
        }
    }
    return meshmodel;
}

void L0ByFaceVector::chooseTypeOfNeighbor()
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

   initFaceSpreadNeighborFace();
}

void L0ByFaceVector::initInfo()
{
    List *faceneighborfacetail = NULL;
    int i, j;
    vector<int>vec;

    arinfocnt = 0;
    infocnt = (int)meshmodel->numtriangles;
    info = new Info*[infocnt];

    if(NULL == info){
        printf("error at L0ByFaceVector::initInfo-----\n");
        return ;
    }

    for(i = 0;i < infocnt - arinfocnt;i++){
        faceneighborfacetail = faceNeighborFace[i];
        vec.clear();
        while(faceneighborfacetail){
            vec.push_back(faceneighborfacetail->data);
            faceneighborfacetail = faceneighborfacetail->next;
        }
        info[i + arinfocnt] = new Info[(int)vec.size()];
        if(NULL == info[i + arinfocnt]){
            printf("error at L0ByFaceVector::initInfo------> info[%d]\n",i + arinfocnt);
            return;
        }

        for(j = 0;j < (int)vec.size();j++){
            info[i + arinfocnt][j].data = vec[j];
            info[i + arinfocnt][j].cnt = (int)vec.size();
            info[i + arinfocnt][j].w = 0.0;
        }
    }
    faceneighborfacetail = NULL;
}

void L0ByFaceVector::updateInfo()
{
    int i, j;
    List *faceneighborfacetail = NULL;

    for(i = 0;i < arinfocnt;i++){
        for(j = 0;j < info[i][0].cnt;j++){
            info[i][j].w = 0;
        }
    }

    vector<double>distv;
    vector<double>normv;
    vector<double>areav;
    vector<double>wv;
    double w;
    int k;

    double max_norm,max_dist;
    double distw, normw, areaw;
    double wsum;

    for(i = 0;i < (int)meshmodel->numtriangles;i++){
        distv.clear();
        normv.clear();
        areav.clear();
        max_dist = max_norm = -1;

        List *faceneighborfacetail = faceNeighborFace[i];
        while(faceneighborfacetail){
            if(i != faceneighborfacetail->data){
                areaw = getAreaByIndex(faceneighborfacetail->data);
                distw = getDistanceFaceCenterToFaceCenter(i, faceneighborfacetail->data);
                normw =  getSignalByFaceIndex(i, faceneighborfacetail->data);

                max_dist = max(max_dist, distw);
                max_norm = max(max_norm, normw);

                areav.push_back(areaw);
                distv.push_back(distw);
                normv.push_back(normw);
            }
            faceneighborfacetail =  faceneighborfacetail->next;
        }
        max_dist = max_dist > 1e-12 ? max_dist : 1e-12;
        max_norm = max_norm > 1e-12 ? max_norm : 1e-12;

        wsum = 0.0;
        wv.clear();
        for(j = 0;j < (int)areav.size();j++){
            w = areav[j] * exp(- distv[j] * distv[j] / (2 * max_dist * max_dist) - normv[j] * normv[j] / (2 * max_norm * max_norm));
            wsum += w;
            wv.push_back(w);
        }
        wsum = wsum > 1e-12 ? wsum : 1e-12;

        faceneighborfacetail = faceNeighborFace[i];
        k = 0;
        j = 0;
        for(j = 0;j < info[i][0].cnt;j++){
            if(i != faceneighborfacetail->data){
                info[i][j].w = wv[k] / wsum;
//                info[i][j].w = 1.0 / (int)areav.size();
                k++;
            }else{
                info[i][j].w = -1;
            }
            faceneighborfacetail = faceneighborfacetail->next;
        }
    }
    faceneighborfacetail = NULL;
}

void L0ByFaceVector::delInfo()
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

void L0ByFaceVector::delFaceNeighborFace()
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

void L0ByFaceVector::initFaceSpreadNeighborFace()
{
    if(faceNeighborFace == NULL){
        printf("error in L0ByFaceVector::initFaceSpreadNeighborFace--->\n");
        return;
    }
    List *tail1,*tail2;
    int i;
    vector<int>v;
    int num;
    faceSpreadNeighborFace = new List*[(int)meshmodel->numtriangles];
    for(i = 0;i < (int)meshmodel->numtriangles;i++){
        faceSpreadNeighborFace[i] = NULL;
    }

    for(i = 0;i < (int)meshmodel->numtriangles;i++){
        tail1 = faceNeighborFace[i];
        v.clear();

        while(tail1){
            tail2 = faceNeighborFace[tail1->data];
            while(tail2){
                v.push_back(tail2->data);
                tail2 = tail2->next;
            }
            tail1 = tail1->next;
        }

        std::sort(v.begin(), v.end());
        num  = std::unique(v.begin(), v.end()) - v.begin();
        insertList(faceSpreadNeighborFace, i, v, num);
    }
    tail1 = tail2 = NULL;
}

void L0ByFaceVector::delFaceSpreadNeighborFace()
{
    if(faceSpreadNeighborFace)
    {
        for(int i = 0;i < (int)meshmodel->numtriangles;i++){
            if(NULL != faceSpreadNeighborFace[i]){
                removeList(faceSpreadNeighborFace[i]);
            }
            faceSpreadNeighborFace[i] = NULL;
        }
        delete []faceSpreadNeighborFace;
        faceSpreadNeighborFace = NULL;
    }
}

void L0ByFaceVector::recoveryVerticesByFaceVector()
{
/*    RecoveryByL2 myRecovery(meshmodel, verticesvindices, verticestindices);
    myRecovery.changeModel(model);
    myRecovery.changeFaceVector(faceVector);
    myRecovery.doL2(2);

    return ;
*/
    int i, j ,k;
    double *newmodel = new double[(int)meshmodel->numvertices * 3];
    if(NULL == newmodel){
        printf("error in L0ByFaceVector::recoveryVerticesByFaceVector--------\n");
        return ;
    }

    for(k = 0;k < (int)meshmodel->numvertices;k++)
    {
        for(j = 0;j < 3;j++)
        {
            newmodel[3 * k + j] = 0;
        }
    }

    int ftemp = 30;
    int num, index;
    double f[3], c[3], temp;
    IndexList *verticestindicestail = NULL;
    while(ftemp--)
    {
        for(i = 0;i < (int)meshmodel->numvertices;i++)
        {
            num = 0;
            for(j = 0;j < 3;j++)
            {
                f[j] = 0.0;
            }
            verticestindicestail = verticestindices[i + 1];
            while(verticestindicestail)
            {
                index = verticestindicestail->index;
                for(j = 0;j < 3;j++)
                {
                    c[j] = 0.0;
                }
                temp = 0.0;

                for(j = 0;j < 3;j++)
                {
                    for(k = 0;k < 3;k++)
                    {
                        c[k] += model[3 * (meshmodel->triangles[index].vindices[j] - 1) + k] / 3;
                    }
                }
                for(j = 0;j < 3;j++)
                {
                    temp += faceVector[j][index] * (c[j] - model[3 * i + j]);
                }
                for(j = 0;j < 3;j++)
                {
                    f[j] += temp * faceVector[j][index];
                }
                num++;
                verticestindicestail = verticestindicestail->next;
            }
            for(j = 0;j < 3;j++)
            {
                newmodel[3 * i + j] = model[3 * i + j] + f[j] / num;
            }
        }
        for(i = 0;i < (int)meshmodel->numvertices;i++)
        {
            for(j = 0;j < 3;j++)
            {
                model[3 * i + j] = newmodel[3 * i + j];
            }
        }
    }
    if(NULL != newmodel){
        delete []newmodel;
    }
    newmodel = NULL;
}
