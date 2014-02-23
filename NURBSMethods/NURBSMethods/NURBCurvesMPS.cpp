#include "NURBCurvesMPS.h"
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


double Q[1220]; // rand number generator stuff
int indx;
double cc;
double c; /* current CSWB */
double zc;	/* current SWB `borrow` */
double zx;	/* SWB seed1 */
double zy;	/* SWB seed2 */
size_t qlen;/* length of Q array */


void NURBCurvesMPS::InitiateRandNumGenerator (unsigned long x)
{

	
	assert(sizeof (double) >= 8);
	cc = 1.0 / 9007199254740992.0; // inverse of 2^53rd power
	int i;
	size_t qlen = indx = sizeof Q / sizeof Q[0];
	for (i = 0; i < qlen; i++)
		Q[i] = 0;
	double c = 0.0, zc = 0.0,	/* current CSWB and SWB `borrow` */
		zx = 5212886298506819.0 / 9007199254740992.0,	/* SWB seed1 */
		zy = 2020898595989513.0 / 9007199254740992.0;	/* SWB seed2 */
	int j;
	double s, t;	 /* Choose 32 bits for x, 32 for y */
	unsigned long /*x = 123456789,*/ y = 362436069; /* default seeds * /
												/* Next, seed each Q[i], one bit at a time, */

	if(x==0){x = 123456789;}

	for (i = 0; i < qlen; i++)
	{ /* using 9th bit from Cong+Xorshift */
		s = 0.0;
		t = 1.0;
		for (j = 0; j < 52; j++)
		{
			t = 0.5 * t; /* make t=.5/2^j */
			x = 69069 * x + 123;
			y ^= (y << 13);
			y ^= (y >> 17);
			y ^= (y << 5);
			if (((x + y) >> 23) & 1)
				s = s + t; /* change bit of s, maybe */
		}	 /* end j loop */
		Q[i] = s;
	} /* end i seed loop, Now generate 10^9 RandNumGenerator's: */

}
double NURBCurvesMPS::RandNumGenerator()
{ /* Takes 14 nanosecs, Intel Q6600,2.40GHz */
	int i, j;
	double t; /* t: first temp, then next CSWB value */
	/* First get zy as next lag-2 SWB */
	t = zx - zy - zc;
	zx = zy;
	if (t < 0)
	{
		zy = t + 1.0;
		zc = cc;
	}
	else
	{
		zy = t;
		zc = 0.0;
	}
	
	/* Then get t as the next lag-1220 CSWB value */
	if (indx < 1220)
		t = Q[indx++];
	else
	{ /* refill Q[n] via Q[n-1220]-Q[n-1190]-c, */
		for (i = 0; i < 1220; i++)
		{
			j = (i < 30) ? i + 1190 : i - 30;
			t = Q[j] - Q[i] + c; /* Get next CSWB element */
			if (t > 0)
			{
				t = t - cc;
				c = cc;
			}
			else
			{
				t = t - cc + 1.0;
				c = 0.0;
			}
			Q[i] = t;
		}	 /* end i loop */
		indx = 1;
		t = Q[0]; /* set indx, exit 'else' with t=Q[0] */
	} /* end else segment; return t-zy mod 1 */
	return ((t < zy) ? 1.0 + (t - zy) : t - zy);
} /* end RandNumGenerator() */

void NURBCurvesMPS::NURBSDartThrowing(size_t _num_active,double*_active,double _s,double*_tmp_active,  //active cell/segemnts stuff
	                                  size_t&_num_samples,double**&_samples,double _r_input,double _tol, // samples stuff
	                                  double** ctrl_p,double _u_max,size_t *knot,size_t K,size_t N,bool kts)// nurb curve stuff
{
	//srand ( time(NULL) ); // activate for different experiments
	InitiateRandNumGenerator(rand());

	ConstructNURBS nurb;
	size_t num_darts,idart,rand_index,iactive,V,num_tmp_active;
	double RF(0.8),rand_u,uu,xx,yy,u_end,u_st,x_st,y_st,x_end,y_end;


	for (int iref=0;iref<30;iref++){ 

		num_darts=RF*_num_active;		

		//dart throwing
		for(idart=0;idart<num_darts;idart++){
			rand_index=size_t((_num_active-1) * RandNumGenerator()); //gives the u start of active intervel
			rand_u=RandNumGenerator()*_s; //gives rand position on the active intervel

			uu=_active[rand_index]+rand_u; 
			
			nurb.PointOnNURBCurve(uu,N,K,ctrl_p,knot,xx,yy,kts);						
						 
			if( !Conflicting(xx,yy,uu,_num_samples,_samples,_r_input,ctrl_p,_u_max,knot,K,N,kts,_tol) ){
				//add sample if not conflicting
				_samples[_num_samples][0]=xx;
				_samples[_num_samples][1]=yy;
				_samples[_num_samples][2]=uu;
				_num_samples++;
			}
		}

		num_tmp_active=0;
		
		if(false){
			nurb.PlotNURB_ps(K,N,ctrl_p,knot,kts,_u_max);
			for(V=0;V<_num_samples;V++){
				PlotThatPoint(_samples[V][0],_samples[V][1],1,ctrl_p,N,_r_input);
			}
		}

		for(iactive=0;iactive<_num_active;iactive++){

			u_st=_active[iactive];
			u_end=u_st+_s;
			if(u_end>=_u_max){u_end=_u_max-_tol;}

			nurb.PointOnNURBCurve(u_st ,N,K,ctrl_p,knot,x_st ,y_st ,kts);						
			nurb.PointOnNURBCurve(u_end,N,K,ctrl_p,knot,x_end,y_end,kts);

			if(false){
				PlotThatPoint(x_st ,y_st,0,ctrl_p,N,_r_input );
				PlotThatPoint(x_end,y_end,0,ctrl_p,N,_r_input);
			}

			if(Covered(x_st,y_st,u_st,x_end,y_end,u_end,_samples,_num_samples,_r_input,ctrl_p,_u_max,knot,K,N,kts,_tol)){continue;}

			_tmp_active[num_tmp_active]=u_st;
			num_tmp_active++;
			_tmp_active[num_tmp_active]=u_st+(_s/2.0);
			num_tmp_active++;
		}
		if(num_tmp_active==0){break;}
		else{
			_num_active=num_tmp_active;
			for(iactive=0;iactive<num_tmp_active;iactive++){
				_active[iactive]=_tmp_active[iactive];
			}
		}
		_s/=2.0;
	}
}
bool NURBCurvesMPS::Covered(double xx1, double yy1,double uu1,double xx2, double yy2,double uu2,double**_samples,size_t _num_samples,double _r_input,
	                       double** ctrl_p,double _u_max,size_t *knot,size_t K,size_t N,bool kts,double _tol)
{	
	//TODO use kd-tree
	size_t V;
	double dist1,dist2;
	for(V=0;V<_num_samples;V++){
		dist1=DistOnNurb(_samples[V][2],uu1,ctrl_p,_u_max,knot,K,N,kts,_tol);
		dist2=DistOnNurb(_samples[V][2],uu2,ctrl_p,_u_max,knot,K,N,kts,_tol);
		if(dist1<_r_input && dist2<_r_input){return true;}
	}
	return false;
}
double NURBCurvesMPS::DistOnNurb(double u1, double u2,double** ctrl_p,double _u_max,size_t *knot,size_t K,size_t N,bool kts,double _tol)
{
	//divide the distance into n interval, and measure euclidean distance between 
	//each two consecutive points
	double umin,umax,u_s,x_prv,y_prv,xx,yy,u,dist(0);
	size_t n(100),V;
	if(u1<u2){umax=u2;umin=u1;}
	else{umax=u1;umin=u2;}
	u_s=double(umax-umin)/double(n);
	bool open;

	if(u2==u1){return 0.0;}

	if(Dist(ctrl_p[1][0],ctrl_p[1][1],0,ctrl_p[N][0],ctrl_p[N][1],0)<_tol*_tol){open=false;}
	else{open=true;}
	
	ConstructNURBS nurb;

	nurb.PointOnNURBCurve(umin,N,K,ctrl_p,knot,x_prv,y_prv,kts);	
	for(V=1;V<=n;V++){
		u=umin+V*u_s;
		if(u>=_u_max){u=_u_max-_tol;}
		nurb.PointOnNURBCurve(u,N,K,ctrl_p,knot,xx,yy,kts);
		dist+=Dist(xx,yy,0,x_prv,y_prv,0);
		x_prv=xx;
		y_prv=yy;
	}

	//if it's open curve, we have just one bath
	//if it's closed curve, then the reverse bath might gives shorter distance
	if(!open){
		bool part1=false;//
		double dist1(0),x_nxt,y_nxt,uu,u_nxt;

		nurb.PointOnNURBCurve(umax,N,K,ctrl_p,knot,xx,yy,kts);
		uu=umax;
		while(true){
			u_nxt=uu+u_s;

			if(u_nxt>=_u_max){
				uu=0;
				nurb.PointOnNURBCurve(uu,N,K,ctrl_p,knot,xx,yy,kts);
				part1=true;
				continue;
			}
			if(u_nxt>umin && part1){break;}
				
			nurb.PointOnNURBCurve(u_nxt,N,K,ctrl_p,knot,x_nxt,y_nxt,kts);
			dist1+=Dist(xx,yy,0,x_nxt,y_nxt,0);

			xx=x_nxt;
			yy=y_nxt;
			uu=u_nxt;
		}


		if(dist1<dist){return dist1;}
	}

	
	return dist;


}
double NURBCurvesMPS::Dist(double x1,double y1, double z1, double x2, double y2, double z2)
{
	double dx,dy,dz;
	dx=x1-x2;
	dy=y1-y2;
	dz=z1-z2;
	dx*=dx;
	dy*=dy;
	dz*=dz;

	return sqrt(dx+dy+dz);
}
void NURBCurvesMPS::PlotThatPoint(double xx, double yy,size_t disk,double** ctrl_p,size_t N,double _r_input)
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
bool NURBCurvesMPS::Conflicting(double xx, double yy,double uu,size_t _num_samples,double**_samples,double _r_input,double**ctrl_p,double _u_max,size_t*knot,size_t K, size_t N,bool kts,double _tol)
{	
	//TODO use kd-tree
	size_t V;
	double dist;
	for(V=0;V<_num_samples;V++){
		dist=DistOnNurb(_samples[V][2],uu,ctrl_p,_u_max,knot,K,N,kts,_tol);
		if(dist<_r_input){return true;}
	}
	return false;
}