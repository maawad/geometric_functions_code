// Geomertical functions code
// Authors: Muhammad A. Awad & Ahmed H. Mahmoud
// Started On: 2-12-2014
// 

#include <iostream>
#include "GeometricalFunctions.h"

using namespace std;

GeometricalFunctions mylib;

void test_d_spheres_intersections()
{
	size_t num_dim = 4;
	mylib._spheres = new double*[num_dim];

	for (size_t i = 0; i < num_dim; i++)
		mylib._spheres[i] = new double[num_dim+1];

	size_t* spheres = new size_t[num_dim];
	double* q1 = new double[num_dim];
	double* q2 = new double[num_dim];

	for (size_t i = 0; i < num_dim; i++)
		spheres[i] = i;
	
	// sphere 0
	mylib._spheres[0][0] = 0.0;
	mylib._spheres[0][1] = 0.3;
	mylib._spheres[0][2] = 0.6;
	mylib._spheres[0][3] = 0.4;
	mylib._spheres[0][4] = 1.0;

	// sphere 1
	mylib._spheres[1][0] = 0.9;
	mylib._spheres[1][1] = 1.1;
	mylib._spheres[1][2] = 1.4;
	mylib._spheres[1][3] = 0.9;
	mylib._spheres[1][4] = 1.0;

	// sphere 2	
	mylib._spheres[2][0] = 0.0;
	mylib._spheres[2][1] = 0.9;
	mylib._spheres[2][2] = 0.9;
	mylib._spheres[2][3] = 1.2;
	mylib._spheres[2][4] = 1.0;

	// sphere 3
	mylib._spheres[3][0] = 0.9;
	mylib._spheres[3][1] = 0.3;
	mylib._spheres[3][2] = 0.9;
	mylib._spheres[3][3] = 1.3;
	mylib._spheres[3][4] = 1.0;

	///////////////// TEST /////////////

	mylib.get_d_spheres_intersections(num_dim,spheres,q1,q2);

	////////////////////////////////////

	// clear memory
	for (size_t i = 0; i < num_dim; i++)
		delete [] mylib._spheres[i];
	delete [] mylib._spheres;

	delete [] spheres;
	delete [] q1;
	delete [] q2;

}

void main()
{
	printf("Welcome !! \n");

	// test functions goes here

	test_d_spheres_intersections();


	system("pause");
}


