#include "stdafx.h"
#include "Cov.h"

#include "MathFunctions.h"
#include <math.h>
#include <time.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include <iostream>

#include <fstream>

using namespace std;

Cov::Cov()
{
	length = 9;
	n = 0;
	weightd = 1.0;
	weightn = 1.0;
	weightl = 1.0;

	verOfNeighborVer = NULL;
	trianglesOfNeighborVer = NULL;

	facetnorm = NULL;
	vernormal = NULL;
	laplaciancoordinates = NULL;

	featureVector = NULL;

	vectorOfCoV = NULL;

	weightCov = NULL;
	delta = NULL;
}


Cov::~Cov()
{
	if (featureVector)
	{
		delete []featureVector;
		featureVector = NULL;
	}

	if (facetnorm)
	{
		delete []facetnorm;
		facetnorm = NULL;
	}

	if (vernormal)
	{
		delete []vernormal;
		vernormal = NULL;
	}

	if (laplaciancoordinates)
	{
		delete []laplaciancoordinates;
		laplaciancoordinates = NULL;
	}

	if (vectorOfCoV)
	{
		delete []vectorOfCoV;
		vectorOfCoV = NULL;
	}

	if (delta)
	{
		delete []delta;
		delta = NULL;
	}

	Deletes(verOfNeighborVer, trianglesOfNeighborVer, weightCov, n);
}

void Cov::ReturnNumofVerofMeshodel(GLMmodel *meshmodel)
{
	n = (int)meshmodel->numvertices;
}


void Cov::FindNeighborVerTriangles(GLMmodel *meshmodel)
{
	if (!verOfNeighborVer)
	{
		verOfNeighborVer = new IndexList * [meshmodel->numvertices+1];
		verOfNeighborVer[0] = NULL;
	}

	if (!trianglesOfNeighborVer)
	{
		trianglesOfNeighborVer = new IndexList * [meshmodel->numvertices+1];
		trianglesOfNeighborVer[0] = NULL;
	}


	for (int i = 1; i <= (int)meshmodel->numvertices; i++)
	{
		verOfNeighborVer[i] = new IndexList;
		verOfNeighborVer[i]->index = 0;
		verOfNeighborVer[i]->next = NULL;

		trianglesOfNeighborVer[i] = new IndexList;
		trianglesOfNeighborVer[i]->index = 0;
		trianglesOfNeighborVer[i]->next = NULL;

	}

	IndexList *p1, pIndexList;
	int indexv, indexvo, m[2];

	for (int j = 0; j < (int)meshmodel->numtriangles; j++)
	{
		for (int k = 0; k < 3; k++)
		{
			indexv = meshmodel->triangles[j].vindices[k];

			p1 = new IndexList;
			p1->index = j;
			p1->next = trianglesOfNeighborVer[indexv]->next;
			trianglesOfNeighborVer[indexv]->next = p1;

			trianglesOfNeighborVer[indexv]->index++;

			if (k == 0)
			{
				m[0] = 1;
				m[1] = 2;
			}
			else if (k == 1)
			{
				m[0] = 0;
				m[1] = 2;
			}
			else
			{
				m[0] = 0;
				m[1] = 1;
			}

			for (int l = 0; l < 2; l++)
			{
				indexvo = meshmodel->triangles[j].vindices[m[l]];

				if (!pIndexList.IsInList(indexvo, verOfNeighborVer[indexv]))  //判断是否有相同的邻点
				{
					p1 = new IndexList;
					p1->index = indexvo;
					p1->next = verOfNeighborVer[indexv]->next;
					verOfNeighborVer[indexv]->next = p1;

					verOfNeighborVer[indexv]->index++;

				}
			}
		}
	}
}


void Cov::ComputeFaceNormal(GLMmodel *meshmodel)
{
	if (!facetnorm)
	{
		facetnorm = new double[3*(meshmodel->numtriangles+1)];
	}

	int indexv1, indexv2, indexv3;
	double uvec[3], vvec[3];

	for (int i = 0; i < (int)meshmodel->numtriangles; i++)
	{
		indexv1 = meshmodel->triangles[i].vindices[0];
		indexv2 = meshmodel->triangles[i].vindices[1];
		indexv3 = meshmodel->triangles[i].vindices[2];

		for (int j = 0; j < 3; j++)
		{
			uvec[j] = meshmodel->vertices[3*indexv2+j] - meshmodel->vertices[3*indexv1+j];
			vvec[j] = meshmodel->vertices[3*indexv3+j] - meshmodel->vertices[3*indexv1+j];
		}

		CrossProd(uvec, vvec, facetnorm+3*(i+1));  //叉积
		Normalize(facetnorm+3*(i+1));

		meshmodel->triangles[i].findex = i+1;
	}
}


void Cov::ComputeVertexNormal(GLMmodel *meshmodel)
{
	if (!vernormal)
	{
		vernormal = new double[3*(meshmodel->numvertices+1)];
	}

	for (int i = 1; i <= (int)meshmodel->numvertices; i++)
	{
		double sum[3] = {0, 0, 0};

		IndexList *p = trianglesOfNeighborVer[i]->next;

		while (p)
		{
			LinearOperation(sum, facetnorm+3*meshmodel->triangles[p->index].findex, 1.0);
			p = p->next;
		}

		for (int l = 0; l < 3; l++)
		{
			vernormal[3*i + l] = sum[l]/trianglesOfNeighborVer[i]->index;
		}

		Normalize(vernormal+3*i);

	}
}


bool Cov::IsInTriangle(GLMmodel *pmesh, int indexv1, int indexv2, 
	                       int *indexv3, int indext)
{
	bool flag = true;

	if (((int)pmesh->triangles[indext].vindices[0] == indexv1 && (int)pmesh->triangles[indext].vindices[1] == indexv2)
		||((int)pmesh->triangles[indext].vindices[0] == indexv2 && (int)pmesh->triangles[indext].vindices[1] == indexv1))
	{
		*indexv3 = (int)pmesh->triangles[indext].vindices[2];
	}
	else if (((int)pmesh->triangles[indext].vindices[0] == indexv1 && (int)pmesh->triangles[indext].vindices[2] == indexv2)
		||((int)pmesh->triangles[indext].vindices[0] == indexv2 && (int)pmesh->triangles[indext].vindices[2] == indexv1))
	{
		*indexv3 = (int)pmesh->triangles[indext].vindices[1];
	}
	else if (((int)pmesh->triangles[indext].vindices[1] == indexv1 && (int)pmesh->triangles[indext].vindices[2] == indexv2)
		||((int)pmesh->triangles[indext].vindices[1] == indexv2 && (int)pmesh->triangles[indext].vindices[2] == indexv1))
	{
		*indexv3 = (int)pmesh->triangles[indext].vindices[0];
	}
	else
	{
		flag = false;
	}

	return flag;
}

void Cov::ComputeCotangentWeight(double *pcotangentweight, int indexv, 
	int *indexcw, GLMmodel *pmesh, IndexList **pverOfNeighborVer, IndexList **ptrianglesOfNeighborVer)
{
	int indexlr, k;
	float dot, crossprod[3], sum = 0.0;
	IndexList *p2, *p3;

	p2 = pverOfNeighborVer[indexv]->next;
	while (p2)
	{
		pcotangentweight[*indexcw] = 0.0;
		k = 0;

		p3 = ptrianglesOfNeighborVer[indexv]->next;
		while (p3)
		{
			if (IsInTriangle(pmesh, indexv, p2->index, &indexlr, p3->index))
			{
				dot = Dot(pmesh->vertices+3*indexlr, pmesh->vertices+3*indexv, pmesh->vertices+3*p2->index);
				CrossProd(pmesh->vertices+3*indexlr, pmesh->vertices+3*indexv, pmesh->vertices+3*p2->index, crossprod);

				pcotangentweight[*indexcw] += dot/Length(crossprod);

				k++;

				if (k == 2)
				{
					break;
				}
			}

			p3 = p3->next;
		}

		sum += pcotangentweight[*indexcw];

		p2 = p2->next;
		(*indexcw)++;
	}

	for (k = (*indexcw-pverOfNeighborVer[indexv]->index); k < *indexcw; k++)
	{
		pcotangentweight[k] /= sum;
	}
}


void Cov::ConstructLaplacianCoordinatesCotangentWeight(double *plc, double *pcw, 
	int *indexcw, int indexv, IndexList **pverOfNeighborVer, float *vertices)
{
	for (int i = 0; i < 3; i++)
	{
		plc[i] = vertices[3*indexv+i];
	}

	IndexList *p1 = pverOfNeighborVer[indexv]->next;
	while (p1)
	{
		for (int i = 0; i < 3; i++)
		{
			plc[i] -= pcw[*indexcw] * vertices[3*p1->index+i];
		}

		p1 = p1->next;
		(*indexcw)++;
	}
}


void Cov::ComputeLaplaciancoordinates(GLMmodel *pmesh, IndexList **pverOfNeighborVer, IndexList **ptrianglesOfNeighborVer)
{
	int numcotangentweight = 0, indexcw;
	double *cotangentweight;
	int i;

	for (i = 1; i <= (int)pmesh->numvertices; i++)
	{
		numcotangentweight += pverOfNeighborVer[i]->index;
	}

	cotangentweight = new double [numcotangentweight];

	indexcw = 0;
	for (i = 1; i <= (int)pmesh->numvertices; i++)
	{
		ComputeCotangentWeight(cotangentweight, i, &indexcw, pmesh, pverOfNeighborVer, ptrianglesOfNeighborVer);
	}

	if (!laplaciancoordinates)
	{
		laplaciancoordinates = new double [3*(pmesh->numvertices+1)];
	}
	
	indexcw = 0;
	for (i = 1; i <= (int)pmesh->numvertices; i++)
	{
		ConstructLaplacianCoordinatesCotangentWeight(laplaciancoordinates+3*i, cotangentweight,
			&indexcw, i, pverOfNeighborVer, pmesh->vertices);
	}

	delete []cotangentweight;
	cotangentweight = NULL;
}


void Cov::GetFeatureVector(GLMmodel *meshmodel)
{
	if (!featureVector)
	{
		featureVector = new double[length*(meshmodel->numvertices+1)];
	}

	FindNeighborVerTriangles(meshmodel);
	ComputeFaceNormal(meshmodel);
	ComputeVertexNormal(meshmodel);
	ComputeLaplaciancoordinates(meshmodel, verOfNeighborVer, trianglesOfNeighborVer);

	double sum1, sum2;

	for (int i = 1; i <= (int)meshmodel->numvertices; i++)
	{
		sum1 = 0;
		sum2 = 0;

		for (int k= 0; k < 3; k++)
		{
			sum1 += meshmodel->vertices[3*i+k]*meshmodel->vertices[3*i+k];
			sum2 += laplaciancoordinates[3*i+k]*laplaciancoordinates[3*i+k];
		}

		for (int l = 0; l < 3; l++)
		{
			if (sum1 < 10e-10)
			{
				sum1 = 10e-10;
			}
			featureVector[length*i+l] = weightd*meshmodel->vertices[3*i+l]/sqrt(sum1);
			
			featureVector[length*i+l+3] = weightn*vernormal[3*i+l];
			
			if (sum2 < 10e-10)
			{
				sum2 = 10e-10;
			}
			featureVector[length*i+l+6] = weightl*laplaciancoordinates[3*i+l]/sqrt(sum2);
		}
	}

}


void Cov::AddVertex(GLMmodel *meshmodel, IndexList **L)
{
	for (int i = 1; i <= (int)meshmodel->numvertices; i++)
	{
		IndexList *p;
		p = L[i]->next;

		while (p)
		{
			IndexList *h;
			h = L[p->index]->next;
			while (h)
			{
				int indexvo;
				indexvo = h->index;
				IndexList *p1, pIndexList;

				if (!pIndexList.IsInList(indexvo, L[i]))
				{
					p1 = new IndexList;
					p1->index = indexvo;
					p1->next = L[i]->next;
					L[i]->next = p1;

					L[i]->index++;
				}

				h = h->next;
			}
			p = p->next;
		}
	}
	
}

void Cov::ComputeNrneighbor(GLMmodel *meshmodel, IndexList **H, int r)
{
	
	if (r > 1)
	{
		for (int i = 1; i < r; i++)
		{
			AddVertex(meshmodel, H);
		}
	}
}


void Cov::ComputeMeanOfFeature(IndexList **neighbor, double *mean, int indexv)
{
	IndexList *p = neighbor[indexv]->next;
	int k ;

	while (p)
	{
		for (k = 0; k < length; k++)
		{
		    mean[k] += featureVector[length*p->index+k]/(1.0 * neighbor[indexv]->index);
		}

		p = p->next;
	}

}

void Cov::ComputeCov(IndexList **neighbor, gsl_matrix *cov, double *mean, int indexv)
{
	double sum;

	for (int i = 0; i < length; i++)
	{
		for (int j = i; j < length; j++)
		{
			sum = 0;
			IndexList *p = neighbor[indexv]->next;
			while (p)
			{
				sum += (featureVector[length*p->index+i]-mean[i]) * (featureVector[length*p->index+j]-mean[j]);

				p = p->next;
			}

			double x;
			x = sum/(neighbor[indexv]->index-1.0);
			gsl_matrix_set(cov, i, j, x);
			gsl_matrix_set(cov, j, i, x);
		}
	}
}


void Cov::Cholesky(gsl_matrix * A)
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

void Cov::GetVectorOfCoV(GLMmodel *meshmodel)
{
	if (!vectorOfCoV)
	{
		vectorOfCoV = new double[(2*length+1)*length*meshmodel->numvertices];
	}

	gsl_matrix *Cov;
	Cov = gsl_matrix_alloc(length, length);

	double *mean = new double[length];
	int i;

	IndexList ** verOfNeighborVer2;
	verOfNeighborVer2 = new IndexList * [meshmodel->numvertices+1];
	verOfNeighborVer2[0] = NULL;
	
	for (i = 1; i <= (int)meshmodel->numvertices; i++)
	{
		verOfNeighborVer2[i] = new IndexList;
		verOfNeighborVer2[i]->index = 0;
		verOfNeighborVer2[i]->next = NULL;
	}

	IndexList *p, *p1, *p2;
	for (i = 1; i <= (int)meshmodel->numvertices; i++)
	{
		p = verOfNeighborVer[i]->next;
		p2 = verOfNeighborVer2[i];

		while (p)
		{
			p1 = new IndexList;
			p1->index = p->index;
			p1->next = NULL;
			p2->next = p1;
			p2 = p2->next;
			
			verOfNeighborVer2[i]->index++;

			p = p->next;
		}
	}

	clock_t start,finish;
	double totaltime;
	
	start = clock();
	ComputeNrneighbor(meshmodel, verOfNeighborVer2, 2);
	finish = clock();
	
	totaltime=(double)(finish-start)/CLOCKS_PER_SEC;
	cout<<totaltime<<endl;

	for (i = 1; i <= (int)meshmodel->numvertices; i++)
	{
		gsl_matrix_set_zero(Cov);
		for (int l = 0; l < length; l++)
		{
			mean[l] = 0;
		}

		ComputeMeanOfFeature(verOfNeighborVer2, mean, i);
		ComputeCov(verOfNeighborVer2, Cov, mean, i);
		Cholesky(Cov);

		for (int i = 0; i < length; i++)
		{
			for (int j = i+1; j < length; j++)
			{
				gsl_matrix_set(Cov, i, j, 0);
			}
		}

		for (int k = 0; k < length; k++)
		{
			vectorOfCoV[(i-1)*(2*length+1)*length+k] = mean[k];
			for (int l = 0; l < length; l++)
			{
				vectorOfCoV[(i-1)*(2*length+1)*length+length*(l+1)+k] = 
					sqrt(2.0)*sqrt(length+0.0)*gsl_matrix_get(Cov,k,l);

				vectorOfCoV[(i-1)*(2*length+1)*length+length*(length+l+1)+k] = 
					-sqrt(2.0)*sqrt(length+0.0)*gsl_matrix_get(Cov,k,l);
			}

		}
	}

	delete []mean;
	mean = NULL;
	gsl_matrix_free(Cov);
	Cov = NULL;

	if (verOfNeighborVer2)
	{
		for (int i = 1; i <= (int)meshmodel->numvertices; i++)
		{
			p1 = verOfNeighborVer2[i];

			while (p1)
			{
				p2 = p1->next;
				delete p1;
				p1 = p2;
			}

			verOfNeighborVer2[i] = NULL;
		}

		delete []verOfNeighborVer2;
		verOfNeighborVer2 = NULL;
	}
}


void Cov::Computedelta(GLMmodel *meshmodel)
{
	if (!delta)
	{
		delta = new double[meshmodel->numvertices+1];
	}
	double maxd, diff;

	for (int i = 1; i <= (int)meshmodel->numvertices; i++)
	{
		IndexList *p;
		p = verOfNeighborVer[i]->next;

		maxd = 0;
		while (p)
		{
			diff = 0;
			for (int k = 0; k < (2*length+1)*length; k++)
			{
				diff += (vectorOfCoV[(2*length+1)*length*(i-1)+k] - vectorOfCoV[(2*length+1)*length*(p->index-1)+k])*
					(vectorOfCoV[(2*length+1)*length*(i-1)+k] - vectorOfCoV[(2*length+1)*length*(p->index-1)+k]);
			}

			if (maxd < diff)
			{
				maxd = diff;
			}

			p = p->next;
		}

		delta[i] = maxd;
		if (delta[i] < 10e-10)
		{
			delta[i] = 10e-10;
		}
	}

}

void Cov::ComputeWeightCov(GLMmodel *meshmodel)
{
	if (!weightCov)
	{
		weightCov = new double *[meshmodel->numvertices+1];
		weightCov[0] = NULL;

		for (int i = 1; i <= (int)meshmodel->numvertices; i++)
		{
			weightCov[i] = new double[verOfNeighborVer[i]->index];
		}
	}

	double diff;
	for (int i = 1; i <= (int)meshmodel->numvertices; i++)
	{
		IndexList *p;
		p = verOfNeighborVer[i]->next;

		int j = 0;

		while (p)
		{
			diff = 0;
			for (int k = 0; k < (2*length+1)*length; k++)
			{
				diff += (vectorOfCoV[(2*length+1)*length*(i-1)+k] - vectorOfCoV[(2*length+1)*length*(p->index-1)+k])*
					(vectorOfCoV[(2*length+1)*length*(i-1)+k] - vectorOfCoV[(2*length+1)*length*(p->index-1)+k]);
			}

			weightCov[i][j++] = exp(-diff/(2*delta[i]));

			p = p->next;
		}
	}
}

void Cov::Denoising(GLMmodel *meshmodel)
{
	GetFeatureVector(meshmodel);
		cout << "yes" << endl;

	GetVectorOfCoV(meshmodel);
		cout << "yes" << endl;

	Computedelta(meshmodel);
		cout << "yes" << endl;

	ComputeWeightCov(meshmodel);
		cout << "yes" << endl;

	for (int i = 1; i <= (int)meshmodel->numvertices; i++)
	{
		double sum = 0;
		for (int j = 0; j < verOfNeighborVer[i]->index; j++)
		{
			sum += weightCov[i][j];
		}
		if (sum < 10e-10)
		{
			sum = 10e-10;
		}

		double sumx = 0, sumy = 0, sumz = 0;
		IndexList *p;
		p = verOfNeighborVer[i]->next;
		j = 0;
		while (p)
		{
			sumx += weightCov[i][j]*meshmodel->vertices[3*p->index+0]/sum;
			sumy += weightCov[i][j]*meshmodel->vertices[3*p->index+1]/sum;
			sumz += weightCov[i][j]*meshmodel->vertices[3*p->index+2]/sum;

			j++;
			p = p->next;
		}

		meshmodel->vertices[3*i+0] = sumx;
		meshmodel->vertices[3*i+1] = sumy;
		meshmodel->vertices[3*i+2] = sumz;

	}
}

void Cov::Deletes(IndexList **pverOfNeighborVer, IndexList **ptrianglesOfNeighborVer, double **pweightCov, int pn)
{
	IndexList *p1, *p2;

	if (pverOfNeighborVer)
	{
		for (int i = 1; i <= (int)pn; i++)
		{
			p1 = pverOfNeighborVer[i];

			while (p1)
			{
				p2 = p1->next;
				delete p1;
				p1 = p2;
			}

			pverOfNeighborVer[i] = NULL;
		}

		delete []pverOfNeighborVer;
		pverOfNeighborVer = NULL;
	}

	if (ptrianglesOfNeighborVer)
	{
		for (int i = 1; i <= (int)pn; i++)
		{
			p1 = ptrianglesOfNeighborVer[i];

			while (p1)
			{
				p2 = p1->next;
				delete p1;
				p1 = p2;
			}

			ptrianglesOfNeighborVer[i] = NULL;
		}

		delete []ptrianglesOfNeighborVer;
		ptrianglesOfNeighborVer = NULL;
	}

	if (pweightCov)
	{
		for (int i = 1; i <= n; i++)
		{
			delete []pweightCov[i];
		}

		delete []pweightCov;
		pweightCov = NULL;
	}
}