#include "L0.h"
#define debug

L0::L0(GLMmodel *originmeshmodel,IndexList **originverticesvindices,IndexList **originverticestindices)
{
        meshmodel = originmeshmodel;
        verticesvindices = originverticesvindices;
        verticestindices = originverticestindices;

        verNeighborVer = NULL;
        verSpreadNeighborVer = NULL;
        faceFirstNeighborFace = NULL;
        faceSecondNeighborFace = NULL;
        edge = NULL;
        faceVector = NULL;
        verVector = NULL;
        faceAllNeighborFace = NULL;

        initModel();
}

L0::~L0()
{
    if(model){
        delete []model;
        model = NULL;
    }

    if(faceVector){
        delVector(faceVector);
        faceVector = NULL;
    }
    if(verVector){
        delVector(verVector);
        verVector = NULL;
    }

    delEdge();
    delVerNeighborVer();
    delVerSpreadNeighborVer();
    delFaceFirstNeighborFace();
    delFaceSecondNeighborFace();
    delFaceAllNeighborFace();
}

void L0::initModel()
{
    int i, j;
    model = new double[(int)meshmodel->numvertices * 3];

    for(i = 1;i <= (int)meshmodel->numvertices; i++){
        for(j = 0;j < 3;j++)
            model[3 * (i - 1) + j] = meshmodel->vertices[3 * i + j];
    }
}

//一环领域
void L0::initVerNeighborVer()
{
    IndexList *verticesvindicestail = NULL;
    int num;

    verNeighborVer = new List*[(int)meshmodel->numvertices];
    int i;
    if(NULL == verNeighborVer){
        printf("Error in L0::initVerNeighborVer\n");
        return ;
    }
    for(i = 0;i < (int)meshmodel->numvertices;i++){
        verNeighborVer[i] = NULL;
    }
    for(i = 0; i< (int)meshmodel->numvertices;i++){
        num = 0;
        verticesvindicestail = verticesvindices[i + 1];
        while(verticesvindicestail){
            num++;
            verticesvindicestail = verticesvindicestail->next;
        }
        num++;

        int *a = new int[num];
        if(NULL == a){
            printf("Error in L0::initVerNeighbor -- a\n");
            return ;
        }

        num = 0;
        verticesvindicestail = verticesvindices[i + 1];
        while(verticesvindicestail){
            a[num++] = verticesvindicestail->index - 1;
            verticesvindicestail = verticesvindicestail->next;
        }
        a[num++] = i;
        std::sort(a,a+num);
        insertList(verNeighborVer, i, a, num);
    }
}

void L0::delVerNeighborVer()
{
    int i;
	if(verNeighborVer)
	{
		for(i = 0;i < (int)meshmodel->numvertices;i++){
			if(verNeighborVer[i])
				removeList(verNeighborVer[i]);
			verNeighborVer[i] = NULL;
		}
	}
    if(verNeighborVer){
        delete []verNeighborVer;
        verNeighborVer = NULL;
    }
}

//二环领域
void L0::initVerSpreadNeighborVer()
{
    IndexList *verticesvindicestail1 = NULL;
    IndexList *verticesvindicestail2 = NULL;
    int num;
    int *a;
    int i, j;

    verSpreadNeighborVer = new List*[(int)meshmodel->numvertices];
    if(NULL == verSpreadNeighborVer){
        printf("Error in L0::initEdgeNeighborVer\n");
        return ;
    }
    for(i = 0;i < (int)meshmodel->numvertices;i++){
        verSpreadNeighborVer[i] = NULL;
    }
    for(i = 0; i< (int)meshmodel->numvertices;i++){
        num = 0;
        verticesvindicestail1 = verticesvindices[i + 1];
        while(verticesvindicestail1){
            verticesvindicestail2 = verticesvindices[verticesvindicestail1->index];
            while(verticesvindicestail2){
                num++;
                verticesvindicestail2 = verticesvindicestail2->next;
            }
            verticesvindicestail1 = verticesvindicestail1->next;
        }
        if(0 == num){
            printf("error in function L0::initEdgeNeighborver ---- num\n");
            return ;
        }

        a = new int[num];
        j = 0;
        if(NULL == a){
            printf("Error in L0::initVerNeighbor -- a\n");
            return ;
        }

        verticesvindicestail1 = verticesvindices[i + 1];
        while(verticesvindicestail1){
            verticesvindicestail2 = verticesvindices[verticesvindicestail1->index];
            while(verticesvindicestail2)
            {
                a[j++] = verticesvindicestail2->index - 1;
                verticesvindicestail2 = verticesvindicestail2->next;
            }
            verticesvindicestail1 = verticesvindicestail1->next;
        }
        std::sort(a,a + num);
        num = std::unique(a, a + num) - a;
        insertList(verSpreadNeighborVer, i, a, num);
    }
}

void L0::delVerSpreadNeighborVer()
{
    int i;
	if(verSpreadNeighborVer){
		for(i = 0;i < (int)meshmodel->numvertices;i++){
			if(verSpreadNeighborVer[i])
				removeList(verSpreadNeighborVer[i]);
			verSpreadNeighborVer[i] = NULL;
		}
	}
    if(verSpreadNeighborVer){
        delete []verSpreadNeighborVer;
        verSpreadNeighborVer = NULL;
    }
}

//面的一类领域面
void L0::initFaceFirstNeighborFace()
{
    if(NULL != faceFirstNeighborFace){
        delFaceFirstNeighborFace();
    }
    int i, num, j;
    IndexList *verticestindicestail = NULL;
    std::vector<int> v;
    int verindex;

    faceFirstNeighborFace = new List*[(int)meshmodel->numtriangles];
    if(NULL == faceFirstNeighborFace){
        printf("error at L0::initFaceFirstNeighborFace------->\n");
        return ;
    }
    for(i = 0;i < (int)meshmodel->numtriangles;i++){
        faceFirstNeighborFace[i] = NULL;
    }

    for(i = 0;i < (int)meshmodel->numtriangles;i++){
        v.clear();

        for(j = 0; j < 3;j++){
            verindex = meshmodel->triangles[i].vindices[j];
            verticestindicestail = verticestindices[verindex];

            while(verticestindicestail){
                if(1 == checkCommonVerticeOfNeighborFace(i, verticestindicestail->index)){
                    v.push_back(verticestindicestail->index);
                }
                verticestindicestail = verticestindicestail->next;
            }
        }

        v.push_back(i);
        std::sort(v.begin(), v.end());
        num  = std::unique(v.begin(), v.end()) - v.begin();
        insertList(faceFirstNeighborFace, i, v, num);
    }
}

void L0::delFaceFirstNeighborFace()
{
    if(faceFirstNeighborFace)
    {
        for(int i = 0;i < (int)meshmodel->numtriangles;i++){
            if(NULL != faceFirstNeighborFace[i]){
                removeList(faceFirstNeighborFace[i]);
            }
            faceFirstNeighborFace[i] = NULL;
        }
        delete []faceFirstNeighborFace;
        faceFirstNeighborFace = NULL;
    }
}

//面的二类领域面
void L0::initFaceSecondNeighborFace()
{
    if(NULL != faceSecondNeighborFace){
        delFaceSecondNeighborFace();
    }
    int i, num, j;
    IndexList *verticestindicestail = NULL;
    std::vector<int> v;
    int verindex;

    faceSecondNeighborFace = new List*[(int)meshmodel->numtriangles];
    if(NULL == faceSecondNeighborFace){
        printf("error at L0::initFaceSecondNeighborFace-->\n");
        return ;
    }
    for(i = 0;i < (int)meshmodel->numtriangles;i++){
        faceSecondNeighborFace[i] = NULL;
    }

    for(i = 0;i < (int)meshmodel->numtriangles;i++){
        v.clear();

        for(j = 0; j < 3;j++){
            verindex = meshmodel->triangles[i].vindices[j];
            verticestindicestail = verticestindices[verindex];

            while(verticestindicestail){
                if(2 == checkCommonVerticeOfNeighborFace(i, verticestindicestail->index)){
                    v.push_back(verticestindicestail->index);
                }
                verticestindicestail = verticestindicestail->next;
            }
        }

        v.push_back(i);
        std::sort(v.begin(), v.end());
        num  = std::unique(v.begin(), v.end()) - v.begin();
        insertList(faceSecondNeighborFace, i, v, num);
    }
}

void L0::delFaceSecondNeighborFace()
{
    if(faceSecondNeighborFace)
    {
        for(int i = 0;i < (int)meshmodel->numtriangles;i++){
            if(NULL != faceSecondNeighborFace[i]){
                removeList(faceSecondNeighborFace[i]);
            }
            faceSecondNeighborFace[i] = NULL;
        }
        delete []faceSecondNeighborFace;
        faceSecondNeighborFace = NULL;
    }
}

//面的所有领域面
void L0::initFaceAllNeighborFace()
{
    int i, num, j;
    IndexList *verticestindicestail = NULL;
    std::vector<int> v;
    int verindex;

    if(NULL != faceAllNeighborFace){
        delFaceAllNeighborFace();
    }
    faceAllNeighborFace = new List*[(int)meshmodel->numtriangles];
    if(NULL == faceAllNeighborFace){
        printf("error at initFaceAllNeighborFace------>\n");
        return ;
    }
    for(i = 0;i < (int)meshmodel->numtriangles;i++){
        faceAllNeighborFace[i] = NULL;
    }

    for(i = 0;i < (int)meshmodel->numtriangles;i++){
        v.clear();

        for(j = 0; j < 3;j++){
            verindex = meshmodel->triangles[i].vindices[j];
            verticestindicestail = verticestindices[verindex];

            while(verticestindicestail){
                if(1 == checkCommonVerticeOfNeighborFace(i, verticestindicestail->index)){
                    v.push_back(verticestindicestail->index);
                }
                else if(2 == checkCommonVerticeOfNeighborFace(i, verticestindicestail->index)){
                    v.push_back(verticestindicestail->index);
                }
                verticestindicestail = verticestindicestail->next;
            }
        }

        v.push_back(i);
        std::sort(v.begin(), v.end());
        num  = std::unique(v.begin(), v.end()) - v.begin();
        insertList(faceAllNeighborFace, i, v, num);
    }
    verticestindicestail = NULL;
}

void L0::delFaceAllNeighborFace()
{
    if(faceAllNeighborFace)
    {
        for(int i = 0;i < (int)meshmodel->numtriangles;i++){
            if(NULL != faceAllNeighborFace[i]){
                removeList(faceAllNeighborFace[i]);
            }
            faceAllNeighborFace[i] = NULL;
        }
        delete []faceAllNeighborFace;
        faceAllNeighborFace = NULL;
    }
}



void L0::initFaceVector()
{
    int i;
    if(faceVector){
        delVector(faceVector);
        faceVector = NULL;
    }

    faceVector = new double *[3];
    if(NULL == faceVector){
        printf("error in L0::initFaceVector\n");
        return ;
    }

    for(i = 0;i < 3;i++){
        faceVector[i] = new double [(int)meshmodel->numtriangles];
        if(NULL == faceVector[i]){
            printf("error in L0::initFaceVector ------->faceVector[%d]\n",i);
            return ;
        }
    }
}

void L0::initVerVector()
{
    int i;
    if(verVector){
        delVector(verVector);
        verVector = NULL;
    }
    verVector = new double *[3];
    if(NULL == verVector){
        printf("error in L0::initFaceVector\n");
        return ;
    }
    for(i = 0;i < 3;i++){
        verVector[i] = new double [(int)meshmodel->numtriangles];
        if(NULL == verVector[i]){
            printf("error in L0::initFaceVector ------->verVector[%d]\n",i);
            return ;
        }
    }
}

void L0::updateFaceVector()
{
    int tindex,i ;
    double vec[3];

    for(tindex = 0;tindex < (int)meshmodel->numtriangles;tindex++)
    {
        getFaceVectorByFaceIndex(vec, tindex);
        for(i = 0;i < 3;i++){
            faceVector[i][tindex] = vec[i];
        }
    }
}

void L0::updateVerVectorByArea()
{
    bool faceVectorFlag = false;
    int verindex, i;
    double faceArea;
    double m;

    IndexList *ftail = NULL;

    //若没有面法相
    if(faceVector){
        faceVectorFlag = true;

        initFaceVector();
        updateFaceVector();
    }

    for(verindex = 0; verindex < (int)meshmodel->numvertices; verindex++){
        ftail = verticesvindices[verindex];
        for(i = 0;i < 3;i++){
            verVector[i][verindex] = 0;
        }
        m = 0;

        while(ftail){
            faceArea = getAreaByIndex(ftail->index);
            for(i = 0; i < 3; i++){
                verVector[i][verindex] += faceVector[i][ftail->index] * faceArea;
            }
            ftail = ftail->next;
        }

        for(i = 0;i < 3;i++){
            m += verVector[i][verindex] * verVector[i][verindex] ;
        }
        m = sqrt(m);
        if(m < 1e-12){
            m = 1e-12;
        }
        for(i = 0;i < 3;i++){
            verVector[i][verindex] /= m;
        }
    }

    ftail = NULL;
    if(faceVectorFlag){
        delVector(faceVector);
        faceVector = NULL;
    }
}

void L0::delVector(double **pvector)
{
    int i;
    if(pvector){
        for(i = 0;i < 3;i++){
            delete []pvector[i];
            pvector[i] = NULL;
        }
        delete []pvector;
        pvector = NULL;
    }
}

int L0::checkCommonEdge(int x,int y,int commonindex)
{
    int i,j;
    for(i = 0;i < 3;i++){
        if((int)meshmodel->triangles[x].vindices[i] != commonindex){
            for(j = 0;j < 3;j++){
                if((int)meshmodel->triangles[y].vindices[j] != commonindex && meshmodel->triangles[x].vindices[i] == meshmodel->triangles[y].vindices[j]){
                    if(commonindex < (int)meshmodel->triangles[x].vindices[i])
                        return meshmodel->triangles[x].vindices[i];
                }
            }
        }
    }
    return -1;
}

int L0::checkCommonVerticeOfNeighborFace(int afaceindex,int bfaceindex)
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

double L0::getPointDistance(int x, int y)
{
    double pp = 0;
    for(int i = 0;i < 3;i++){
        pp += (model[3 * x + i] - model[3 * y + i])
            * (model[3 * x + i] - model[3 * y + i]);
    }
    return sqrt(pp);
}

double L0::getPointAngle(double z,double x,double y)
{
    double temp = 2 * x * y;
    if(temp < 1e-14)
        temp = 1e-14;

    double cosangle = ((x * x + y * y) - z * z) / ( temp );

    return acos(cosangle);
}

double L0::getPointFactor(int index1,int index2,int index3,int index4)
{
    double a[3],b[3],factor;
    int i ;
    factor = 0.0;
    for(i = 0;i < 3;i++)
    {
        a[i] = model[3 * index1 + i] - model[3 * index2 + i];
        b[i] = model[3 * index3 + i] - model[3 * index4 + i];
    }
    for(i = 0;i < 3;i++)
        factor += a[i] * b[i];
    return factor;
}

double L0::getAreaByIndex(int findex)
{
    double u[3], v[3], r[3];

    int a = meshmodel->triangles[findex].vindices[0] - 1;
    int b = meshmodel->triangles[findex].vindices[1] - 1;
    int c = meshmodel->triangles[findex].vindices[2] - 1;

    u[0] = model[3 * b + 0] - model[3 * a + 0];
    u[1] = model[3 * b + 1] - model[3 * a + 1];
    u[2] = model[3 * b + 2] - model[3 * a + 2];

    v[0] = model[3 * c + 0] - model[3 * a + 0];
    v[1] = model[3 * c + 1] - model[3 * a + 1];
    v[2] = model[3 * c + 2] - model[3 * a + 2];

    r[0] = u[1] * v[2] - u[2] * v[1];
    r[1] = u[2] * v[0] - u[0] * v[2];
    r[2] = u[0] * v[1] - u[1] * v[0];

    return sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
}

double L0::getDistanceFaceCenterToFaceCenter(int afindex,int bfindex)
{
    double a[3], b[3], m;
    int i, j;

    for(i = 0;i < 3;i++){
        a[i] = b[i] = 0;
    }
    for(i = 0;i < 3;i++){
        for(j = 0;j < 3;j++){
            a[j] += model[(meshmodel->triangles[afindex].vindices[i] - 1) * 3 + j] / 3;
            b[j] += model[(meshmodel->triangles[bfindex].vindices[i] - 1) * 3 + j] / 3;
        }
    }

    m = 0.0;
    for(i = 0;i < 3;i++){
        m = (a[i] - b[i]) * (a[i] - b[i]);
    }
    m = sqrt(m);
    m = m > 1e-12 ? m : 1e-12;
    return m;
}

void L0::getEdge()
{
    int i,j;
    IndexList *verticestindicestail11 = NULL;
    IndexList *verticestindicestail22 = NULL;
    Edge *edgetail = NULL;
    edgeNum = 0;

    for(i = 1;i <= (int)meshmodel->numvertices;i++){
        IndexList *verticestindicestail11 = verticestindices[i];
        while(NULL != verticestindicestail11){
            IndexList *verticestindicestail22 = verticestindicestail11->next;
            while(NULL != verticestindicestail22){
                int resultindex = checkCommonEdge((int)verticestindicestail11->index,(int)verticestindicestail22->index,i);
                if(resultindex != -1){
                    if(NULL == edgetail){
                        edge = createEdgeNode();
                        edgetail = edge;
                    }else{
                        edgetail->next = createEdgeNode();
                        edgetail = edgetail->next;
                    }

                    edgetail->p[0] = i - 1;
                    edgetail->p[2] = resultindex - 1;
                    edgeNum++;

                    GLMtriangle trip2index = meshmodel->triangles[verticestindicestail11->index];
                    GLMtriangle trip4index = meshmodel->triangles[verticestindicestail22->index];

                    for(j = 0;j < 3;j++){
                        if((int)trip2index.vindices[j] != i && (int)trip2index.vindices[j] != resultindex){
                            edgetail->p[1] = trip2index.vindices[j] - 1;
                        }
                    }
                    for(j = 0;j < 3;j++){
                        if((int)trip4index.vindices[j] != i && (int)trip4index.vindices[j] != resultindex){
                            edgetail->p[3] = trip4index.vindices[j] - 1;
                        }
                    }
					edgetail->f[0] = min(verticestindicestail11->index,verticestindicestail22->index);
					edgetail->f[1] = max(verticestindicestail11->index,verticestindicestail22->index);
                }
                verticestindicestail22 = verticestindicestail22->next;
            }
            verticestindicestail11 = verticestindicestail11->next;
        }
    }
    verticestindicestail22 = NULL;
    verticestindicestail11 = NULL;
}

void L0::delEdge(){
    if(edge){
        removeList(edge);
        edge = NULL;
    }
}

Edge *L0::createEdgeNode()
{
    Edge *s = new Edge();
    if(s == NULL){
        printf("error at L0::createEdgeNode\n");
    }
    for(int i = 0;i < 4;i++){
        s->p[i] = 0;
    }
    s->next = NULL;
    return s;
}

List *L0::createListNode()
{
    List *s = new List();
    if(s == NULL){
        printf("error at L0::createListNode\n");
    }
    s->data = 0.0;
    s->next = NULL;

    return s;
}

void L0::insertList(List **relation, int place, int *a,int num)
{
    List *tail = relation[place];
    for(int i = 0;i < num;i++){
        if(tail == NULL){
            relation[place] = createListNode();
            tail = relation[place];
        }else{
            tail->next = createListNode();
            tail = tail->next;
        }
        tail->data = a[i];
    }
    tail = NULL;
}
void L0::insertList(List **relation, int place, std::vector<int> v,int num)
{
    List *tail = relation[place];
    for(int iter = 0;iter < num;iter++){
        if(tail == NULL){
            relation[place] = createListNode();
            tail = relation[place];
        }else{
            tail->next = createListNode();
            tail = tail->next;
        }
        tail->data = v[iter];
    }
    tail = NULL;
}

double L0::getAngleByFaceIndex(int afindex,int bfindex)
{
    double avec[3];
    double bvec[3];
    getFaceVectorByFaceIndex(avec, afindex);
    getFaceVectorByFaceIndex(bvec, bfindex);

    double alen = 0.0;
    double blen = 0.0;
    double sum = 0.0;
    int i;

    for(i = 0;i < 3;i++){
        alen += avec[i] * avec[i];
        blen += bvec[i] * bvec[i];
        sum += avec[i] * bvec[i];
    }
    alen = sqrt(alen);
    blen = sqrt(blen);

    alen = alen > 1e-12 ? alen : 1e-12;
    blen = blen > 1e-12 ? blen : 1e-12;

    double r = sum / (alen * blen);
    if(r > 1.0){
        r = 1.0;
    }
    if(r < -1.0){
        r = -1.0;
    }
    return acos(r);
}

double L0::getSignalByFaceIndex(int afindex,int bfindex)
{
    double avec[3];
    double bvec[3];
    getFaceVectorByFaceIndex(avec, afindex);
    getFaceVectorByFaceIndex(bvec, bfindex);

    double sum = 0.0;
    int i;

    for(i = 0;i < 3;i++){
        sum += (avec[i] - bvec[i]) * (avec[i] - bvec[i]);
    }
    sum = sqrt(sum);
    sum = sum > 1e-12 ? sum : 1e-12;
    return sum;
}

void L0::getFaceVectorByFaceIndex(double *vec, int findex)
{
    double u[3], v[3];
    int i;

    int a = meshmodel->triangles[findex].vindices[0] - 1;
    int b = meshmodel->triangles[findex].vindices[1] - 1;
    int c = meshmodel->triangles[findex].vindices[2] - 1;

    u[0] = model[3 * b + 0] - model[3 * a + 0];
    u[1] = model[3 * b + 1] - model[3 * a + 1];
    u[2] = model[3 * b + 2] - model[3 * a + 2];

    v[0] = model[3 * c + 0] - model[3 * a + 0];
    v[1] = model[3 * c + 1] - model[3 * a + 1];
    v[2] = model[3 * c + 2] - model[3 * a + 2];

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
