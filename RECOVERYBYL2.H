#ifndef RECOVERYBYL2_H
#define RECOVERYBYL2_H
#include "L0.h"

class RecoveryByL2 : public L0{
    public:
        Info **info;
        List **verThreeSeparateRelation;
        double *p;
        double *v;

    public:
        RecoveryByL2(GLMmodel *originmeshmodel,IndexList **originverticesvindices,IndexList **originverticestindices);
        ~RecoveryByL2();

        void doL2(double pbeta);

        void changeModel(double *pModel);
        void changeFaceVector(double **pFaceVector);

        void initInfo();
        void updateInfo();
        void delInfo();

        void initVerThreeSeparateRelation();
        void delVerThreeSeparateRelation();
};
#endif
