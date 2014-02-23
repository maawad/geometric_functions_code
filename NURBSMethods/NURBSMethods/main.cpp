#include "ConstructNURBS.h"
#include "DrawSpheres.h"
#include "NURBCurvesMPS.h"
#include "ProjectPointToNURBCurve.h"

#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>  
#include <fstream>
#include <sstream>
#include <time.h>
#include <cassert>

using namespace std;


#define PI 3.1415926536
#define TwoPI 6.2831853072




//#define Surface // define this if working on curved surfaces


double*** ctrl_p_ij; //control points coordinates (x,y,z,wt) followed by their weight
double** ctrl_p; //control points coordinates (x,y,wt) followed by their weight
double _u_max,_w_max;
size_t *knot,K,N,L,M,*knot_u,*knot_w;//knot vector 
bool kts,kts_u,kts_w;




double _r_input,**_samples,_s,*_active,*_tmp_active,_tol;
size_t _num_samples,_num_active,_num_expected,_num_tmp_active;

	


inline void PlotThatPoint(double xx, double yy,size_t disk)
{
	double xmin(10E-10), ymin(10E-10),xmax(-10E-10), ymax(-10E-10),scale;

	for(size_t V=1;V<=N;V++){
		if(ctrl_p[V][0]<xmin){xmin=ctrl_p[V][0];}
		if(ctrl_p[V][1]<ymin){ymin=ctrl_p[V][1];}

		if(ctrl_p[V][0]>xmax){xmax=ctrl_p[V][0];}
		if(ctrl_p[V][1]>ymax){ymax=ctrl_p[V][1];}
	}
	double lx(xmax-xmin),ly(ymax-ymin);
	scale=8.0/(lx);

	fstream file("NURB.ps",ios::app);	
	file << "1.0 1.0 0.0 setrgbcolor" << endl;
	if(disk==0){
		file << "0.0 0.0 0.0 setrgbcolor" << endl;
	}
	//file<< "0.01 setlinewidth"<<endl;
	
	file<< (xx)*scale<<" "<< (yy)*scale<<" "<<0.02*scale<<" dot"<<endl;
	if(disk==1){
		file<< (xx)*scale<<" "<< (yy)*scale<<" "<<_r_input*scale<<" pink_disk"<<endl;
	}
	
}
inline double Dist(double x1,double y1, double z1, double x2, double y2, double z2)
{
	double dx,dy,dz;
	dx=x1-x2;
	dy=y1-y2;
	dz=z1-z2;
	dx*=dx;
	dy*=dy;
	dz*=dz;

	return dx+dy+dz;

}


int main()
{
	
	size_t V,d,i,j,tuna;
		
	cout << "Enter the Radius" << endl;
	cin >> _r_input;
	cout << endl;


	//****************Reading Input File Starts Here****************///

	//File Formats For Curves:
	//------------
	// order
	//# of control points (N)
	// x0 y0 wt0
	// x1 y1 wt1
	// ...	
	//xN-1 yN-1 wtN-1 
	//1 for auto knot vector calc, 0 for inout knot vector
	//# of knots
	//knot1 
	//knot2 
	//knot3
	//...
	//..


	//File Formats For Surfaces:
	//------------
	// order_u
	//# of control points (N) in u
	// order_w
	//# of control points (M) in w
	// x0 y0 z0 wt0
	// x1 y1 z1 wt1
	// ...	
	//xNM-1 yNM-1 wtNM-1	
	 
	//1 for auto knot_u vector calc, 0 for inout knot_u vector
	//# of knots_u
	//knot_u1 
	//knot_u2 
	//knot_u3
	//...
	//..

	//1 for auto knot_w vector calc, 0 for inout knot_w vector
	//# of knots_w
	//knot_w1 
	//knot_w2 
	//knot_w3
	//...
	//..

	cout<<"\n Reading Input File"<<endl;
	ifstream inputfile;
	inputfile.open("batman.txt");
	//inputfile.open("circle.txt");
	//inputfile.open("curve.txt");
	//inputfile.open("surface.txt");
	inputfile>>K; //oder or order_u 
	inputfile>>N; // num of control points / u

#ifdef Surface
	inputfile>>L; // order_w
	inputfile>>M; // num of control points w	
	
	ctrl_p_ij=new double**[N*M];	
	for(i=1;i<=N;i++){
		ctrl_p_ij[i]=new double*[M+1];
		for(j=1;j<=M;j++){
			ctrl_p_ij[i][j]=new double [4];
			inputfile>>ctrl_p_ij[i][j][0]; //x
			inputfile>>ctrl_p_ij[i][j][1]; //y
			inputfile>>ctrl_p_ij[i][j][2]; //z
			inputfile>>ctrl_p_ij[i][j][3]; //weight
		}
	}
	size_t num_kts,in_u,in_w;
	inputfile>>in_u;
	if(in_u==1){kts_u=true;}
	else if(in_u==0){kts_u=false;}
	else{
		cout<<"Error(0).. invalid input file"<<endl;
		system("pause");
	}
	if(!kts_u){
		inputfile>>num_kts;
		if(num_kts!=N+K){
			cout<<"Error(1).. invalid input length of knots vector"<<endl;
			system("pause");
		}
		knot_u=new size_t[num_kts+1];
		for(V=1;V<=num_kts;V++){
			inputfile>>knot_u[V];
		}
	}else{
		knot_u=new size_t[N+K+1];
	}
	
	
	inputfile>>in_w;
	if(in_w==1){kts_w=true;}
	else if(in_w==0){kts_w=false;}
	else{
		cout<<"Error(2).. invalid input file"<<endl;
		system("pause");
	}
	if(!kts_u){
		inputfile>>num_kts;
		if(num_kts!=M+L){
			cout<<"Error(3).. invalid input length of knots vector"<<endl;
			system("pause");
		}
		knot_w=new size_t[num_kts+1];
		for(V=1;V<=num_kts;V++){
			inputfile>>knot_w[V];
		}
	}else{
		knot_w=new size_t[M+L+1];
	}
	

	_u_max=N-K+1;
	_w_max=M-L+1;
#else
	ctrl_p=new double*[N+1];
	for(V=1;V<=N;V++){
		ctrl_p[V]=new double[3];
		inputfile>>ctrl_p[V][0]; //x
		inputfile>>ctrl_p[V][1]; //y
		inputfile>>ctrl_p[V][2]; // weight
	}

	size_t num_kts,in;	
	inputfile>>in;
	if(in==1){kts=true;}
	else if(in==0){kts=false;}
	else{
		cout<<"Error(0).. invalid input file"<<endl;
		system("pause");
	}
	
	if(!kts){
		inputfile>>num_kts;
		if(num_kts!=N+K){
			cout<<"Error(1).. invalid input length of knots vector"<<endl;
			system("pause");
		}
		knot=new size_t[num_kts+1];
		for(V=1;V<=num_kts;V++){
			inputfile>>knot[V];
		}
	}else{
		knot=new size_t[N+K+1];
	}

	//********
	//_u_max=4.0; // use this with circle since it's closed. (still don't know how to calculate u_max for closed nurbs)
	_u_max=N-K+1; // use this with open/any nurbs
	//*********

#endif 
	//****************Reading Input File Ends Here****************///
		
	







	ConstructNURBS nurb;
	double x,y,z,u,w;	
	//Draw your nurbs as a start
#ifdef Surface
	nurb.PlotNURB_surface(K,N,L,M,ctrl_p_ij,knot_u,knot_w,kts_u,kts_w,_u_max,_w_max);	
#else
	nurb.PlotNURB_ps(K,N,ctrl_p,knot,kts,_u_max);
#endif

		



	//****************Initilize Sampling Data Structure Starts Here****************///
	cout<<"\n Initialize"<<endl;	
	_tol=10E-10; //tolerance
	_num_samples=0; //total number of samples
	_num_active=50; //number of active segements (similar to active cells)
	_s=_u_max/_num_active; //spacing 
	_num_expected=10E5;
	_active=new double[_num_expected]; //store the active cells/segments (just the starting parameter u of the segment)
	_tmp_active=new double[_num_expected];
	_active[0]=0.0;
	for(V=1;V<_num_active;V++){
		_active[V]=_active[V-1]+_s;
	}

	_samples=new double*[_num_expected];//samples coordinates+parameter
	for(V=0;V<_num_expected;V++){
		_samples[V]=new double[3]; //(x,y,u)
	}
	//add the 1st control point directly (since the samples should represent the curve)
	_samples[_num_samples][0]=ctrl_p[1][0];
	_samples[_num_samples][1]=ctrl_p[1][1];
	_samples[_num_samples][2]=0;
	_num_samples++;
		
	//add the end control point only if it's not identical with the 1st ctrl point (i.e. not a close curve)
	double dx,dy;
	dx=ctrl_p[1][0]-ctrl_p[N][0];
	dy=ctrl_p[1][1]-ctrl_p[N][1];
	dx*=dx;
	dy*=dy;
	if(Dist(ctrl_p[1][0],ctrl_p[1][1],0,ctrl_p[N][0],ctrl_p[N][1],0)>_tol*_tol){
		_samples[_num_samples][0]=ctrl_p[N][0];
		_samples[_num_samples][1]=ctrl_p[N][1];
		_samples[_num_samples][2]=_u_max;
		_num_samples++;
	}
	//****************Initilize Sampling Data Structure Ends Here****************///










	//****************Sampling Starts Here****************///
	cout<<"\n MPS Sampling"<<endl;
	NURBCurvesMPS nurbs_mps;
	nurbs_mps.NURBSDartThrowing(_num_active,_active,_s,_tmp_active,  //active cell/segemnts stuff
	                            _num_samples,_samples,_r_input,_tol, // samples stuff
	                            ctrl_p,_u_max,knot,K,N,kts);// nurb curve stuff
	cout<<"num_points= "<<_num_samples<<endl;
	//Plot the output samples
	for(V=0;V<_num_samples;V++){
		PlotThatPoint(_samples[V][0],_samples[V][1],1);
	}
	//****************Sampling Ends Here****************///
	








	//****************Testing Getting The Min Distance Starts Here****************///
	ProjectPointToNURBCurve min_dist;
	double dd;
	dd=min_dist.MinDist(2.0,0.0,u,ctrl_p,_u_max,knot,K,N,kts);
	//****************Testing Getting The Min Distance Ends Here****************///


	
	return 0;
}
