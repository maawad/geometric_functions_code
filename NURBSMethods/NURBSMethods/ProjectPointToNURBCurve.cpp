#include "ProjectPointToNURBCurve.h"
#include "ConstructNURBS.h"
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>  
#include <fstream>
#include <sstream>
#include <time.h>
#include <cassert>
using namespace std;

double ProjectPointToNURBCurve::MinDist(double xc, double yc,double&uu,double** ctrl_p,double _u_max,size_t *knot,size_t K,size_t N,bool kts)
{
	// TODO apply this algorithem http://www-sop.inria.fr/members/Gang.Xu/English/paper/CAD08_miniDistance.pdf


	//find the min distance between point (xc,ycs) and NURB curve 
	// return distance, uu as the closest point on the curve to (xx,yy)


	//incrementaly moves along the curve, calculate the distance as we move along
	size_t V,d,num(50.0);
	double s(_u_max/double(num)),u,dist_close(10E10),xx,yy,dist;
	ConstructNURBS nurb;
	
	for(V=0;V<num;V++){
				
		u=double(V)*s;

		nurb.PointOnNURBCurve(u,N,K,ctrl_p,knot,xx,yy,kts);

		dist=Dist(xc,yc,0,xx,yy,0);

		if(dist<dist_close){
			uu=u;
			dist_close=dist;
		}
	}

	if(false){
		//plots the closest points and stright line between it and the test point (xc,yc)
		nurb.PointOnNURBCurve(uu,N,K,ctrl_p,knot,xx,yy,kts);

		PlotThatPoint(xc,yc,xx,yy,ctrl_p,N);
	}

	return sqrt(dist_close);
}
double ProjectPointToNURBCurve::Dist(double x1,double y1, double z1, double x2, double y2, double z2)
{
	double dx,dy,dz;
	dx=x1-x2;
	dy=y1-y2;
	dz=z1-z2;
	dx*=dx;
	dy*=dy;
	dz*=dz;

	return (dx+dy+dz);
}
void ProjectPointToNURBCurve::PlotThatPoint(double xx, double yy,double xc,double yc,double** ctrl_p,size_t N)
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
	file << "0.0 0.0 0.0 setrgbcolor" << endl;
	
	//file<< "0.01 setlinewidth"<<endl;
	
	file<< (xx)*scale<<" "<< (yy)*scale<<" "<<0.02*scale<<" dot"<<endl;

	file<<(xx)*scale<<" "<< (yy)*scale<<" moveto"<<endl;
	file<<(xc)*scale<<" "<< (yc)*scale<<" lineto"<<endl;
	file<<"stroke"<<endl;
	
	
}