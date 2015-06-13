#include "ReCov.h"
#include "stdafx.h"

#include "MathFunctions.h"
#include <math.h>
#include <time.h>
#include <iterator>
#include <algorithm>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include <iostream>

#include <fstream>


using namespace std;

ReCov::ReCov(GLMmodel *pmeshmodel,IndexList **pverticesvindices,IndexList **pverticestindices)
{
    meshmodel = pmeshmodel;
    faceNeighborFace = NULL;
    rFaceNeighborFace = NULL;
    dFaceNeighborFace = NULL;
    verticestindices = pverticestindices;
    verticesvindices = pverticesvindices;

    featureVector = NULL;
    faceNormal = NULL;
    vectorOfCoV = NULL;
    weightCov = NULL;
    delta = NULL;

    length = 0;
}

ReCov::~ReCov()
{
    meshmodel = NULL;
    verticestindices = NULL;
    verticesvindices = NULL;


    delFaceNormal();
    delFeatureVector();
}

void ReCov::getFaceNeighborFace()
{
    if(NULL != faceNeighborFace){
        delFaceNeighborFace();
    }
    int i, num, j;
    IndexList *verticestindicestail = NULL;
    std::vector<int> v;
    int verindex;

    faceNeighborFace = new IndexList*[(int)meshmodel->numtriangles];
    if(NULL == faceNeighborFace){
        printf("error at ReCov::initFaceNeighborFace-->\n");
        return ;
    }
    for(i = 0;i < (int)meshmodel->numtriangles;i++){
        faceNeighborFace[i] = NULL;
    }

    for(i = 0;i < (int)meshmodel->numtriangles;i++){
        v.clear();

        for(j = 0; j < 3;j++){
            verindex = meshmodel->triangles[i].vindices[j];
            verticestindicestail = verticestindices[verindex];

            while(verticestindicestail){
                /*if(2 == checkCommonVerticeOfNeighborFace(i, verticestindicestail->index)){
                    v.push_back(verticestindicestail->index);
                }*/
	            if(1 == checkCommonVerticeOfNeighborFace(i, verticestindicestail->index)){
                    v.push_back(verticestindicestail->index);
                }
                verticestindicestail = verticestindicestail->next;
            }
        }

        v.push_back(i);
        std::sort(v.begin(), v.end());
        num  = std::unique(v.begin(), v.end()) - v.begin();
        insertIndexList(faceNeighborFace, i, v, num);
    }
}

void ReCov::delFaceNeighborFace()
{

}

int ReCov::checkCommonVerticeOfNeighborFace(int afaceindex,int bfaceindex)
{
    int i, j;
    if(afaceindex == bfaceindex){
        return 0;
    }
    int a[3], b[3], sum;
    sum = 0;
    for(i = 0;i < 3;i++){
        a[i] = meshmodel->triangles[afaceindex].vindices[i];
        b[i] = meshmodel->triangles[bfaceindex].vindices[i];
    }
    for(i = 0; i < 3;i++){
        for(j = 0;j < 3;j++){
            if(a[i] == b[j]){
                sum++;
            }
        }
    }
    return sum;
}

void ReCov::initDFaceNeighborFace()
{
    if(NULL == dFaceNeighborFace){
        delDFaceNeighborFace();
        dFaceNeighborFace = NULL;
    }
    dFaceNeighborFace = new IndexList*[(int)meshmodel->numtriangles];
    if(NULL == dFaceNeighborFace){
        printf("Errors in ReCov::initRFaceNeighborFace\n");
        return ;
    }
    for(int i = 0;i < (int)meshmodel->numtriangles;i++){
        dFaceNeighborFace[i] = NULL;
    }
}

void ReCov::delDFaceNeighborFace()
{

}

void ReCov::getFaceInsideFaces(IndexList **Head,IndexList **faceNeighborFace,int r)
{
    if(NULL == faceNeighborFace){
        printf("Errors in ReCov::getFaceInsideFaces\n");
        printf("faceNeighborFace is NULL\n");
        return ;
    }

    int i, t;
    std::map<int, char> mp;
    std::map<int, char> mptemp1;
    std::map<int, char> mptemp2;
    std::map<int, char>::iterator it;

    IndexList *tail = NULL;

    for(i = 0;i < (int)meshmodel->numtriangles;i++){
        t = r;
        map<int, char>().swap(mp);
        map<int, char>().swap(mptemp1);
        map<int, char>().swap(mptemp2);

        mp[i] = 1;
        mptemp2[i] = 1;

        while(t--){
            mptemp1 = mptemp2;
            map<int, char>().swap(mptemp2);

            for(it = mptemp1.begin();it != mptemp1.end();it++){
                tail = faceNeighborFace[it->first];
                while(tail){
                    if(!mp[tail->index]){
                        mptemp2[tail->index] = 1;
                    }
                    tail = tail->next;
                }
            }
            for(it = mptemp2.begin();it != mptemp2.end();it++){
                mp[it->first] = 1;
            }
        }
        insertIndexList(Head, i, mp, (int)mp.size());
    }
    tail = NULL;
}

void ReCov::initFaceNormal()
{
    if(NULL != faceNormal){
        delFaceNormal();
        faceNormal = NULL;
    }
    int i;
    faceNormal = new double*[(int)meshmodel->numtriangles];
    if(NULL == faceNormal){
        printf("Errors in ReCov::initFaceNormal\n");
        return ;
    }
    for(i = 0;i < (int)meshmodel->numtriangles;i++){
        faceNormal[i] = new double[3];
        if(NULL == faceNormal[i]){
            printf("Errors in ReCov::initFaceNormal\n");
            return ;
        }
    }
}

void ReCov::updateFaceNormal()
{
    int tindex, i;
    double vec[3];

    for(tindex = 0;tindex < (int)meshmodel->numtriangles;tindex++)
    {
        getFaceVectorByFaceIndex(vec, tindex);
        for(i = 0;i < 3;i++){
            faceNormal[tindex][i] = vec[i];
        }
    }
}

void ReCov::getFaceVectorByFaceIndex(double *vec, int findex)
{
    double u[3], v[3];
    int i;

    int a = meshmodel->triangles[findex].vindices[0];
    int b = meshmodel->triangles[findex].vindices[1];
    int c = meshmodel->triangles[findex].vindices[2];

    u[0] = meshmodel->vertices[3 * b + 0] - meshmodel->vertices[3 * a + 0];
    u[1] = meshmodel->vertices[3 * b + 1] - meshmodel->vertices[3 * a + 1];
    u[2] = meshmodel->vertices[3 * b + 2] - meshmodel->vertices[3 * a + 2];

    v[0] = meshmodel->vertices[3 * c + 0] - meshmodel->vertices[3 * a + 0];
    v[1] = meshmodel->vertices[3 * c + 1] - meshmodel->vertices[3 * a + 1];
    v[2] = meshmodel->vertices[3 * c + 2] - meshmodel->vertices[3 * a + 2];

    vec[0] = u[1] * v[2] - u[2] * v[1];
    vec[1] = u[2] * v[0] - u[0] * v[2];
    vec[2] = u[0] * v[1] - u[1] * v[0];

    double sum = 0.0;
    for(i = 0;i < 3;i++){
        sum += vec[i] * vec[i];
    }
    sum = sqrt(sum);
    sum = sum > 1e-12 ? sum : 1e-12;
    for(i = 0;i < 3;i++){
        vec[i] /= sum;
    }
}

void ReCov::delFaceNormal()
{
    int i;
    if(NULL != faceNormal){
        for(i = 0;i < (int)meshmodel->numtriangles;i++){
            delete[] faceNormal[i];
            faceNormal[i] = NULL;
        }
        delete[] faceNormal;
        faceNormal = NULL;
    }
}

void ReCov::getFeatureVector()
{
    double mid[3];
    length = 6;
    int i, j, k, vidx;
    featureVector = new double[(int)meshmodel->numtriangles * length];
    if(NULL == featureVector){
        printf("Errors in ReCov::getFeatureVector\n");
        return ;
    }

    initFaceNormal();
    updateFaceNormal();

    for(i = 0;i < (int)meshmodel->numtriangles;i++){
        for(j = 0;j < 3;j++){
            featureVector[length * i + j] = faceNormal[i][j];
        }
        for(j = 0;j < 3;j++){
            mid[j] = 0.0;
        }
        for(j = 0;j < 3;j++){
            vidx = meshmodel->triangles[i].vindices[j];
            for(k = 0;k < 3;k++){
                mid[k] += meshmodel->vertices[3 * vidx + k];
            }
        }
        for(j = 0;j < 3;j++){
            mid[j] /= 3;
            featureVector[length * i + j + 3] = mid[j];
        }
    }

    delFaceNormal();
}

void ReCov::delFeatureVector()
{

}

void ReCov::computeMeanOfFeature(IndexList **dFaceNeighborFace,double *mean,int vindex)
{
    int i, num;
    for(i = 0;i < length;i++){
        mean[i] = 0;
    }
    num = 0;

    IndexList *tail = dFaceNeighborFace[vindex];
    while(tail){
        for(i = 0;i < length;i++)
            mean[i] += featureVector[tail->index * length + i];
        num++;
        tail = tail->next;
    }
    for(i = 0;i < length;i++){
        mean[i] /= (1.0 * num);
    }
    tail = NULL;
}

void ReCov::computeCov(IndexList **dFaceNeighborFace,gsl_matrix *Cov,double *mean, int vindex)
{
    double sum;
    int num;
    IndexList *tail = NULL;

    for (int i = 0; i < length; i++)
    {
        for (int j = i; j < length; j++)
        {
            sum = 0;
            num = 0;
            tail = dFaceNeighborFace[vindex];
            while (tail)
            {
                num++;
                sum += (featureVector[length*tail->index+i]-mean[i]) * (featureVector[length*tail->index+j]-mean[j]);
                tail = tail->next;
            }

            double x;
            x = sum/(1.0 * (num - 1));
            gsl_matrix_set(Cov, i, j, x);
            gsl_matrix_set(Cov, j, i, x);
        }
    }
    tail = NULL;
}

void ReCov::getVectorOfCoV()
{
    if (!vectorOfCoV)
    {
        vectorOfCoV = new double[(2*length+1)*length*meshmodel->numtriangles];
    }
    int i, j, k;
    gsl_matrix *Cov;
    Cov = gsl_matrix_alloc(length, length);
    double *mean = new double[length];
    if(NULL == mean){
        printf("Errors in ReCov::getVectorOfCoV\n");
        return ;
    }

    for(i = 0;i < (int)meshmodel->numtriangles;i++){
        gsl_matrix_set_zero(Cov);
        for(j = 0;j < length;j++){
            mean[j] = 0;
        }
        computeMeanOfFeature(dFaceNeighborFace, mean, i);
        computeCov(dFaceNeighborFace, Cov, mean, i);
        cholesky(Cov);
        for (k = 0; k < length;k++)
        {
            for (int o = k+1; o < length; o++)
            {
                gsl_matrix_set(Cov, k, o, 0);
            }
        }

        for (k = 0; k < length; k++)
        {
            vectorOfCoV[(i)*(2*length+1)*length+k] = mean[k];
            for (int l = 0; l < length; l++)
            {
                vectorOfCoV[(i)*(2*length+1)*length+length*(l)+k] =
                    sqrt(2.0)*sqrt(length+0.0)*gsl_matrix_get(Cov,k,l);

                vectorOfCoV[(i)*(2*length+1)*length+length*(length+l)+k] =
                    -sqrt(2.0)*sqrt(length+0.0)*gsl_matrix_get(Cov,k,l);
            }
        }
    }

    delete []mean;
    mean = NULL;
}

void ReCov::delVectorOfCov()
{

}

void ReCov::getDelta()
{
    int i, k;
    IndexList *tail = NULL;
    if (!delta)
    {
        delta = new double[(int)meshmodel->numtriangles];
    }
    double maxd, diff;

    for (i = 0; i < (int)meshmodel->numtriangles; i++)
    {
        tail = dFaceNeighborFace[i];

        maxd = 0;
        while (tail)
        {
            diff = 0;
            for (k = 0; k < (2*length+1)*length; k++)
            {
                diff += (vectorOfCoV[(2*length+1)*length*(i)+k] - vectorOfCoV[(2*length+1)*length*(tail->index)+k])*
                    (vectorOfCoV[(2*length+1)*length*(i)+k] - vectorOfCoV[(2*length+1)*length*(tail->index)+k]);
            }

            if (maxd < diff)
            {
                maxd = diff;
            }

            tail = tail->next;
        }

        delta[i] = maxd;
        if (delta[i] < 10e-10)
        {
            delta[i] = 10e-10;
        }
		delta[i] = 0.2;
    }
    tail = NULL;
}

void ReCov::delDelta()
{

}

void ReCov::getWeightCov()
{
    if(NULL == dFaceNeighborFace){
        printf("Erros in ReCov::getWeightCov\ndFaceNeightborFace is NULL\n");
        return ;
    }
    vector<double> v;
    int i, j, k, point;
    IndexList *tail = NULL;
    double sum;

    if (!weightCov)
    {
        weightCov = new double *[meshmodel->numtriangles];
        for(i = 0;i < meshmodel->numtriangles;i++){
            weightCov[i] = NULL;
        }
    }

    double diff;
    for (i = 0; i < (int)meshmodel->numtriangles; i++)
    {
        vector<double>().swap(v);
        v.clear();
        sum = 0.0;

        tail = dFaceNeighborFace[i];
        j = 0;
        while (tail)
        {
            diff = 0;
            for (k = 0; k < (2*length+1)*length; k++)
            {
                diff += (vectorOfCoV[(2*length+1)*length*(i)+k] - vectorOfCoV[(2*length+1)*length*(tail->index)+k])*
                    (vectorOfCoV[(2*length+1)*length*(i)+k] - vectorOfCoV[(2*length+1)*length*(tail->index)+k]);
            }
            v.push_back(exp(-diff/(2*delta[i])));
            if(i != tail->index){
                sum += exp(-diff/(2*delta[i]));
            }else{
                point = j;
            }
            tail = tail->next;
            j++;
        }
        weightCov[i] = new double[(int)v.size()];

        if(NULL == weightCov[i]){
            printf("errors in getWeightCov\n");
            return ;
        }
        sum = sum > 1e-3 ? sum : 1e-3;

        for(j = 0;j < v.size();j++){
            if(j != point){
                weightCov[i][j] = v[j] / sum;
            }else{
                weightCov[i][j] = v[j];
            }
        }
    }
    tail = NULL;
}

void ReCov::delWeightCov()
{
}

void ReCov::cholesky(gsl_matrix * A)
{
    int M = A->size1;
    int N = A->size2;

    if (M != N)
    {
        cout<<"matrix must be positive definite"<<endl;
    }
    else
    {
        int i,j,k;
        int status = 0;

        /* Do the first 2 rows explicitly.  It is simple, and faster.  And
         * one can return if the matrix has only 1 or 2 rows.
         */

        double A_00 = gsl_matrix_get (A, 0, 0);

        if (A_00 < 10e-10)
        {
            A_00 = 10e-10;
        }

        double L_00 = sqrt(A_00);

        if (A_00 < 0)
        {
            status = GSL_EDOM ;
        }

        gsl_matrix_set (A, 0, 0, L_00);

        if (M > 1)
        {
            double A_10 = gsl_matrix_get (A, 1, 0);
            double A_11 = gsl_matrix_get (A, 1, 1);

            double L_10 = A_10 / L_00;
            double diag = A_11 - L_10 * L_10;

            if (diag < 10e-10)
            {
                diag = 10e-10;
            }

            double L_11 = sqrt(diag);

            if (diag < 0)
            {
                status = GSL_EDOM;
            }

            gsl_matrix_set (A, 1, 0, L_10);
            gsl_matrix_set (A, 1, 1, L_11);
        }

        for (k = 2; k < M; k++)
        {
            double A_kk = gsl_matrix_get (A, k, k);

            for (i = 0; i < k; i++)
            {
                double sum = 0;

                double A_ki = gsl_matrix_get (A, k, i);
                double A_ii = gsl_matrix_get (A, i, i);

                gsl_vector_view ci = gsl_matrix_row (A, i);
                gsl_vector_view ck = gsl_matrix_row (A, k);

                if (i > 0) {
                    gsl_vector_view di = gsl_vector_subvector(&ci.vector, 0, i);
                    gsl_vector_view dk = gsl_vector_subvector(&ck.vector, 0, i);

                    gsl_blas_ddot(&di.vector, &dk.vector, &sum);
                }

                A_ki = (A_ki - sum) / A_ii;
                gsl_matrix_set (A, k, i, A_ki);
            }

            {
                gsl_vector_view ck = gsl_matrix_row (A, k);
                gsl_vector_view dk = gsl_vector_subvector (&ck.vector, 0, k);

                double sum = gsl_blas_dnrm2(&dk.vector);
                double diag = A_kk - sum * sum;

                if (diag < 10e-10)
                {
                    diag = 10e-10;
                }

                double L_kk = sqrt(diag);

                if (diag < 0)
                {
                    status = GSL_EDOM;
                }

                gsl_matrix_set (A, k, k, L_kk);
            }
        }

        /* Now copy the transposed lower triangle to the upper triangle,
         * the diagonal is common.
         */

        for (i = 1; i < M; i++)
        {
            for (j = 0; j < i; j++)
            {
                double A_ij = gsl_matrix_get (A, i, j);
                gsl_matrix_set (A, j, i, A_ij);
            }
        }

        if (status == GSL_EDOM)
        {
            cout<<"matrix must be positive definite"<<endl;
        }
    }
}

void ReCov::insertIndexList(IndexList **H,int place,std::map<int,char> &mp, int num)
{
    IndexList *tail = H[place];
    IndexList *s = NULL;

    std::map<int, char>::iterator it;

    for(it = mp.begin();it != mp.end();it++){
        s = new IndexList();
        s->index = it->first;
        s->next = NULL;

        if(tail == NULL){
            tail = H[place] = s;
        }else{
            tail->next = s;
            tail = tail->next;
        }
    }
    tail = NULL;
}

void ReCov::insertIndexList(IndexList **H,int place,std::vector<int> &v, int num)
{
    IndexList *tail = H[place];
    IndexList *s = NULL;

    for(int i = 0;i < num;i++){
        s = new IndexList();
        s->index = v[i];
        s->next = NULL;

        if(tail == NULL){
            tail = H[place] = s;
        }else{
            tail->next = s;
            tail = tail->next;
        }
    }
    tail = NULL;
}

void ReCov::delIndexList(IndexList **Head,int num)
{

}

void ReCov::doCov()
{
    //特征向量
    getFeatureVector();

    //面片r环领域
    getFaceNeighborFace();
    initDFaceNeighborFace();
    getFaceInsideFaces(dFaceNeighborFace, faceNeighborFace, 1);

    //分解得到的向量
    getVectorOfCoV();
    //获得delta
    getDelta();
    //获得权重
    getWeightCov();
}