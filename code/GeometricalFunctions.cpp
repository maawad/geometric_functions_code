#include "GeometricalFunctions.h"
#include <cstdlib>
#include <iostream>

GeometricalFunctions::GeometricalFunctions(void)
{
}


GeometricalFunctions::~GeometricalFunctions(void)
{
}


/*********************
	Random methods
*********************/
double GeometricalFunctions::generate_a_random_number()
{
	return rand()/RAND_MAX;
}


/*********************
	linear system solvers methods
*********************/
bool GeometricalFunctions::LU_decompose_pivot(double** A,int* pivot,size_t N)
{

	// Crout's LU with partial pivoting
	// ref. <http://ezekiel.vancouver.wsu.edu/~cs330/lectures/linear_algebra/LU.pdf>
	
	int n ;

	int i,j,k,jp1,dum;
	double* tmp;
	double sum,p;
/*		
	for( i = 0 ; i < _d ; i++)
	{
		for( j = 0 ; j < _d ; j++)
		{
			_L[i][j] = 0;
			_U[i][j] = 0;


		}
	}  
	*/
	
	for(i = 0 ; i < N ; i++) pivot[i] = i;
	int d = 1;

	for (j = 0; j < N; j++) 
	{
		for(i = 0; i <= j; i++) 
		{
			sum=0;
			for(k = 0; k < i; k++) 
				sum += (A[i][k] * A[k][j]);  

			A[i][j] -= sum;
		}
		p = fabs(A[j][j]);
		n = j;
		jp1 = j + 1;
		for(i = jp1 ; i < N; i++)
		{
			sum=0;
			for(k = 0; k < j; k++) 
				sum += (A[i][k] * A[k][j]);  

			A[i][j] -= sum;
			if(fabs(A[i][j]) > p)
			{
				p = fabs(A[i][j]);
				n = i;
			}
		}

		if(p < 0)	printf("**Error: Singular Matrix! \n");
		if( n != j)
		{
			tmp = A[n];
			A[n] = A[j];
			A[j] = tmp;
			d *= -1;

			dum = pivot[n];
			pivot[n] = pivot[j];
			pivot[j] = dum;
		}
		for(i = j + 1 ; i < N ; i++)
			A[i][j] /= A[j][j]; 
	}

	return false;

/*
	if(false)	// debugging
	{
		for( i = 0 ; i < _d - 1 ; i++)
		{
			for( j = 0 ; j < _d - 1 ; j++)
			{
				cout << L[i][j] <<"	";


			}
			cout << endl;
		}
			cout << endl;
			cout << endl;
			cout << endl;
			cout << endl;
	
		for( i = 0 ; i < _d - 1 ; i++)
		{
			for( j = 0 ; j < _d - 1 ; j++)
			{
				cout << U[i][j] <<"	";


			}
			cout << endl;
		}
	}
	*/

}
void GeometricalFunctions::solver_pivot(double** a,double* b,int* pivot,size_t N,double* x)
{
	double sum = 0;
	int i,j;

	x[0] = b[pivot[0]];

	for(i = 1 ; i < N ; i++)
	{
		sum = 0 ;

		for(j = 0 ; j < i  ; j++)
			sum += x[j] * a[i][j]; 

		x[i] = (b[pivot[i]] - sum);
	}

	x[N-1] /= a[N-1][N-1];

	for(i = int(N)-2 ; i >= 0 ; i--)
	{
		sum = 0 ;

		for(j = i + 1 ; j < N  ; j++)
			sum += x[j] * a[i][j]; 

		x[i] = (x[i] - sum)/a[i][i];
	}
}



/*********************
	Geometrical funcs
*********************/
bool GeometricalFunctions::get_d_spheres_intersections(size_t num_dim,size_t* spheres_indices,double* q1,double* q2)
{
	// prepare memory
	_planes_normals = new double*[num_dim];
	_B = new double[num_dim];
	_ortho_basis = new double*[num_dim];
	_ortho_tmp = new double[num_dim];
	_I = new double[num_dim];
	_p1 = new double[num_dim]; // point of the d-planes intersectio

	size_t ip0 (spheres_indices[0]),ip1;

	for (size_t isphere = 0; isphere < num_dim - 1; isphere++)
	{
		_planes_normals[isphere] = new double[num_dim];

		ip1 = spheres_indices[isphere + 1];
	
		double dot(0.0);

		for (size_t idim = 0; idim < num_dim; idim++)
		{
			_planes_normals[isphere][idim] = _spheres[ip1][idim] - _spheres[ip0][idim];
			double xm = 0.5 * (_spheres[ip1][idim] + _spheres[ip0][idim]);
			dot += _planes_normals[isphere][idim]*xm;
		}

		_B[isphere] = dot;
	}

	// get the direction of the d vector normal to the d-1 vectors

	// orthonormal basis 
	for (size_t k = 0; k < num_dim - 1; k++)
	{
		_ortho_basis[k] = new double[num_dim];
		
		for (size_t idim = 0; idim < num_dim; idim++)
			_ortho_basis[k][idim] = _planes_normals[k][idim];

		for (size_t j = 0; j < k; j++)
		{
			double dot = 0.0;
			for (size_t idim = 0; idim < num_dim; idim++)
				dot += _ortho_basis[j][idim]*_planes_normals[k][idim];

			_ortho_tmp[j] = dot;

			for (size_t idim = 0; idim < num_dim; idim++)
				_ortho_basis[k][idim] -= _ortho_tmp[j]*_ortho_basis[j][idim];	
		}
		
		double sum = 0.0;
		for (size_t idim = 0; idim < num_dim; idim++)
			sum += _ortho_basis[k][idim]*_ortho_basis[k][idim];

		// normailze
		sum = sqrt(sum);
		for (size_t idim = 0; idim < num_dim; idim++)
			_ortho_basis[k][idim]/=sum;
	}

	
	// random vector as the d vector , then normalize
	
	_planes_normals[num_dim-1] = new double[num_dim];	
	while (true)
	{
		for (size_t idim = 0; idim < num_dim; idim++)
			_planes_normals[num_dim-1][idim] =  generate_a_random_number();

		size_t k = num_dim - 1;
		_ortho_basis[k] = new double[num_dim];
		
		for (size_t idim = 0; idim < num_dim; idim++)
			_ortho_basis[k][idim] = _planes_normals[k][idim];

		for (size_t j = 0; j < k; j++)
		{
			double dot = 0.0;
			for (size_t idim = 0; idim < num_dim; idim++)
				dot += _ortho_basis[j][idim]*_planes_normals[k][idim];

			_ortho_tmp[j] = dot;

			for (size_t idim = 0; idim < num_dim; idim++)
				_ortho_basis[k][idim] -= _ortho_tmp[j]*_ortho_basis[j][idim];	
		}
		
		double sum = 0.0;
		for (size_t idim = 0; idim < num_dim; idim++)
			sum += _ortho_basis[k][idim]*_ortho_basis[k][idim];

		// normailze
		sum = sqrt(sum);
		if(sum < 1E-10) continue; // try again 
		for (size_t idim = 0; idim < num_dim; idim++)
		{
			_ortho_basis[k][idim]/=sum;
			_I[idim] = _ortho_basis[k][idim];
		}
		break; // done
	}

	double* tmp_ptr = _planes_normals[num_dim - 1];
	_planes_normals[num_dim - 1] = _ortho_basis[num_dim - 1];
	_ortho_basis[num_dim - 1] = tmp_ptr;
	
	double dot = 0.0;
	for (size_t idim = 0; idim < num_dim; idim++)
		dot += _spheres[ip0][idim]*_planes_normals[num_dim - 1][idim];
	_B[num_dim - 1] = dot;

	_pivot = new int[num_dim];
	LU_decompose_pivot(_planes_normals,_pivot,num_dim);
	solver_pivot(_planes_normals,_B,_pivot,num_dim,_p1);

	// line sphere intersection 

	// u = -(I.(o-c)) +- sqrt( (I.(o-c))^2 - (o-c)^2 + r^2 )
	// o line first point
	// I line unit vector
	// c sphere centre
	double a(0.0),b(0.0),c(0.0);
	for (size_t idim = 0; idim < num_dim; idim++)
	{
		double fc = _p1[idim] - _spheres[spheres_indices[0]][idim]; //(o - c)
		a += fc * _I[idim];  // (o - c).I
		b += fc*fc;          // (o - c).(o - c)
	}
	c = a*a; // a2
	double r = _spheres[spheres_indices[0]][num_dim];
	double disc = c - b + r*r;

	if(disc < 0.0) // no intersection between the d-spheres
	{
		// clear memory
		delete [] _B;
		delete [] _I;
		delete [] _p1;
		delete [] _ortho_tmp;
		delete [] _pivot;

		for (size_t idim = 0; idim < num_dim; idim++)
		{
			delete [] _planes_normals[idim];
			delete [] _ortho_basis[idim];
		}

		return false;	
	} 

	// intersection point between the d-spheres exist
 	a*= -1;
	double disc_sq = sqrt(disc);
	double u1 = a + disc_sq;
	double u2 = a - disc_sq;

	for (size_t idim = 0; idim < num_dim; idim++)
	{
		q1[idim] = _p1[idim] + u1 * _I[idim];
		q2[idim] = _p1[idim] + u2 * _I[idim];
	}


	// clear memory
	delete [] _B;
	delete [] _I;
	delete [] _p1;
	delete [] _ortho_tmp;
	delete [] _pivot;

	for (size_t idim = 0; idim < num_dim; idim++)
	{
		delete [] _planes_normals[idim];
		delete [] _ortho_basis[idim];
	}
	return true;
}

