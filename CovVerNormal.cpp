#include "CovVerNormal.h"

CovVerNormal::CovVerNormal(GLMmodel *originmeshmodel, IndexList **originverticesvindices, IndexList **originverticestindices) \
        :L0(originmeshmodel, originverticesvindices, originverticestindices)
{
	verWeight = NULL;	
}

CovVerNormal::~CovVerNormal()
{
}