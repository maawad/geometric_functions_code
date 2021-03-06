#include "ConstructSpline.h"
#include "DrawSpheres.h"
//#include "SplineCurvesMPS.h"
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>  
#include <fstream>
#include <sstream>
#include <time.h>
#include <cassert>
using namespace std;

void ConstructSpline::KnotVector(size_t K, size_t N,size_t*&knot)
{
	// open knot vector
		
	size_t i;
	//knot = new size_t[K+N+1];
	for(i=1;i<=K+N;i++){
		if(i<=K && i>=1){knot[i]=0;}
		else if(i>=K+1 && i<=N){knot[i]=i-K;}
		else if (i>=N+1 && i<=N+K){knot[i]=N-K+1;}
		else {
			cout<<"Error at KnotVector()"<<endl;
			system("pause");
		}
		if(i>1){
			if(knot[i-1]>knot[i]){
				cout<<"Error(0) at KnotVector()"<<endl;
				// should be nondecreasing order
				system("pause");
			}
		}
	}
}
size_t ConstructSpline::N_i_1(double u,size_t i,size_t*knot)
{		
	if(u>=knot[i] && u<knot[i+1]){
		return 1;
	}else{
		return 0;
	}
}
double ConstructSpline::BasisFunction(size_t i, size_t K, size_t N,double u, size_t*knot)
{
	if(K<1 || K>N){
		cout<<"Error(2) at BasisFunction()  "<<endl;
		system("pause");
	}
	if(K==1){
		return N_i_1(u,i,knot);
	}
	else{
		double term1,term2,term3,term4,term5,term6;
		term1=(u-knot[i])*BasisFunction(i,K-1,N,u,knot);
		term2=knot[i+K-1]-knot[i];

		term3=(knot[i+K]-u)*BasisFunction(i+1,K-1,N,u,knot);
		term4=knot[i+K]-knot[i+1];

		if(term1==0&&term2==0){term5=0;}
		else{term5=term1/term2;}
		if(term2==0&&term1!=0){
			cout<<"Error(0) at BasisFunction()  "<<endl;
			system("pause");
		}

		if(term3==0&&term4==0){term6=0;}
		else{term6=term3/term4;}
		if(term4==0&&term3!=0){
			cout<<"Error(1) at BasisFunction()  "<<endl;
			system("pause");
		}
		if(term5+term6<0.0){
			cout<<"Error(2) at BasisFunction()  "<<endl;
			system("pause");
		}
		return term5+term6;
	}
}
void ConstructSpline::PointOnSplineCurve(double u,size_t N,size_t K,double** ctrl_p,size_t*&knot,double&x,double&y,bool kts,double _tol,double _u_max)
{
	//http://www.enseignement.polytechnique.fr/informatique/INF555/Slides/lecture4_v2.pdf
	//Open Rational B-Spline
	//0<=u<=N-K+2

	if(kts){KnotVector(K,N,knot);}

	if(u==0.0){u+=_tol;}
	if(u==_u_max){u-=_tol;}

	if(K<=1 || K>N){
		cout<<"Error(0) at PointOnSplineCurve(). Invalid K value"<<endl;
		system("pause");
	}
	

	if(u<0 || u>_u_max){
		cout<<"Error(1) at PointOnSplineCurve(). Invalid u value"<<endl;
		system("pause");
	}

	size_t i,d;
	double bs_fun;

	double term1_x(0),term2_x(0),term1_y(0),term2_y(0),test(0);
	for(i=1;i<=N;i++){
		bs_fun=BasisFunction(i,K,N,u,knot);
		test+=bs_fun;
		if(bs_fun<0){
			cout<<"Error(2) at PointOnSplineCurve()"<<endl;
			system("pause");
		}
		term1_x+=ctrl_p[i][2]*ctrl_p[i][0]*bs_fun;
		term2_x+=ctrl_p[i][2]*bs_fun;		

		term1_y+=ctrl_p[i][2]*ctrl_p[i][1]*bs_fun;
		term2_y+=ctrl_p[i][2]*bs_fun;
	}

	if(abs(test-1.0)>0.00000001 && kts){
		cout<<"Error(3) at PointOnSplineCurve()"<<endl;
		system("pause");
	}
	term2_x=1.0;
	term2_y=1.0;
	x=term1_x/term2_x;
	y=term1_y/term2_y;
}
void ConstructSpline::PlotSpline_ps(size_t K,size_t N,double**ctrl_p,size_t*knot,bool kts,double max_u,double _tol)
{
	fstream file("Spline.ps", ios::out);
	file.precision(30);
	//fstream file("MPS1.ps", ios::app);

	file.precision(30);
	file << "%!PS-Adobe-3.0" << endl;  
	file << "60 60 scale     % one unit = one Centimeter" << endl;

	double xmin(10E-10), ymin(10E-10),xmax(-10E-10), ymax(-10E-10);
    double scale_x, scale_y;
    double shift_x, shift_y;
	size_t V;
	for(V=1;V<=N;V++){
		if(ctrl_p[V][0]<xmin){xmin=ctrl_p[V][0];}
		if(ctrl_p[V][1]<ymin){ymin=ctrl_p[V][1];}

		if(ctrl_p[V][0]>xmax){xmax=ctrl_p[V][0];}
		if(ctrl_p[V][1]>ymax){ymax=ctrl_p[V][1];}
	}
	double lx(xmax-xmin),ly(ymax-ymin);
	

	scale_x=8.0/(lx);
	
	scale_y=10.0/(ly);


	//shift_x = 1.5*scale_x+xmin;	//with circle
	//shift_x = 1.1*scale_x+xmin;	//with open curve 
	shift_x=1.0; //with batman
	shift_y = 0.05 * (28.0 - (1.0) * scale_x)+5.0;	

    file << shift_x << " " << shift_y << " translate" << std::endl;

    file<<"/Times-Roman findfont"<<endl;
	file<<"0.2 scalefont"<<endl;
	file <<"setfont"<<endl;
	
#pragma region Definitions of Shapes:
	file << "/seg      % stack: x1 y1 x2 y2" << std::endl;
    file << "{newpath" << std::endl; 
    file << " moveto" << std::endl;
    file << " lineto" << std::endl;
    file << " closepath" << std::endl;
    file << " 0.01 setlinewidth" << std::endl;
	file << " 0 0 0 setrgbcolor" << std::endl;
    file << " stroke" << std::endl;
    file << "} def" << std::endl;


    file << "/quad      % stack: x1 y1 x2 y2 x3 y3 x4 y4" << std::endl;
    file << "{newpath" << std::endl; 
    file << " moveto" << std::endl;
    file << " lineto" << std::endl;
    file << " lineto" << std::endl;
    file << " lineto" << std::endl;
    file << " closepath" << std::endl;
	//file <<" 0.0 0 0.0 setrgbcolor" <<endl;
    file << " 0.004 setlinewidth" << std::endl;
    file << " stroke" << std::endl;
    file << "} def" << std::endl;


	file <<"/cell  % stack: x1 y1 x2 y2 x3 y3 x4 y4 " <<endl;
	file <<"{" <<endl;
	file <<" newpath" <<endl;
	file <<" moveto" <<endl;
	file <<" lineto" <<endl;
	file <<" lineto" <<endl;
	file <<" lineto" <<endl;
	file <<"closepath" <<endl;
	file <<" 0.5 0.5 0.5 setrgbcolor" <<endl;
	file <<"fill"<<endl;
	file<<"} def" <<endl;


	file << "/quad_bold      % stack: x1 y1_ x2 y2 x3 y3 x4 y4" << endl;
    file << "{newpath" << endl; 
    file << " moveto" << endl;
    file << " lineto" << endl;
    file << " lineto" << endl;
    file << " lineto" << endl;
    file << " closepath" << endl;
    file << " 0.01 setlinewidth" << endl;
    file << " stroke" << endl;
    file << "} def" << endl;


	file << "/pink_disk      % stack: x y r" << endl;
	file << "{newpath" << endl; 
	file << " 0 360 arc closepath" << endl;
	file << "  1.0 0.0 1.0 setrgbcolor" << endl;
	//file << " fill" << endl;
	file << " 0.008 setlinewidth" << endl;
	file << " stroke" << endl;
	file << "} def" << endl;

	file << "/dot      % stack: x y r" << endl;
	file << "{newpath" << endl; 
	file << " 0 360 arc closepath" << endl;	
	file << " fill" << endl;
	file << " 0.01 setlinewidth" << endl;
	file << " stroke" << endl;
	file << "} def" << endl;



   file<<"/quad_white      % stack: x1 y1 x2 y2 x3 y3 x4 y4"<<endl;
   file<<"{newpath"<<endl;
   file<<"moveto"<<endl;
   file<<"lineto"<<endl;
   file<<"lineto"<<endl;
   file<<"lineto"<<endl;
   file<<"closepath"<<endl;
   file<<"gsave"<<endl;
   file<<"1.0 setgray fill"<<endl;
   file<<" grestore"<<endl;
   file<<"} def"<<endl;

#pragma endregion
	
 
   double xx,yy,u;

   
   for(V=1;V<=N;V++){
	   file << "1.0 0.0 0.0 setrgbcolor" << endl;
	    file<< (ctrl_p[V][0])*scale_x<<" "<< (ctrl_p[V][1])*scale_x<<" "<<0.005*scale_x<<" dot"<<endl;
		if(V<N){
			
			//file<<(ctrl_p[V][0])*scale_x<<" "<< (ctrl_p[V][1])*scale_x<<" ";
			//file<<(ctrl_p[V+1][0])*scale_x<<" "<< (ctrl_p[V+1][1])*scale_x<<" seg"<<endl;
			
		}
   }


   file << "0.0 0.0 0.83 setrgbcolor" << endl;
   for(u=0.00000001;u<max_u;u+=0.001){
	   PointOnSplineCurve(u,N,K,ctrl_p,knot,xx,yy,kts,_tol,max_u);
	   
	   file<< (xx)*scale_x<<" "<< (yy)*scale_x<<" "<<0.001*scale_x<<" dot"<<endl;
   }

}

