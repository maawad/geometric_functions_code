#pragma once
class GeometricalFunctions
{
public:
	GeometricalFunctions(void);
	~GeometricalFunctions(void);



/*********************
	Random methods
*********************/

	double generate_a_random_number();


/*********************
	linear system solvers methods
*********************/

	bool LU_decompose_pivot(double** A,int* pivot,size_t N);
	void solver_pivot(double** a,double* b,int* pivot,size_t N,double* x);


/*********************
	Geometrical funcs
*********************/
	bool get_d_spheres_intersections(size_t num_dim,size_t* spheres_indices,double* q1,double* q2);


////////////////////////////////////////////////////////////////
/// Variables
////////////////////////////////////////////////////////////////

			
	// d-balls intersection
	double** _spheres;
	double** _planes_normals; // vectors normal to the planes
	double* _B; // linear system b for solving d planes eqs  
	double** _ortho_basis; // baisis vectors
	double* _ortho_tmp;
	double* _I; // line unit vector
	double* _p1; // intersection point of the d-planes
	double* _q1; // intersection point
	double* _q2; // intersection point
	int* _pivot;
};

