#include "NURBSurfacesMPS.h"
#include "ConstructNURBS.h"
#include "DrawSpheres.h"
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>  
#include <fstream>
#include <sstream>
#include <time.h>
#include <cassert>
using namespace std;


extern	double Q[1220]; // rand number generator stuff
extern	int indx;
extern	double cc;
extern	double c; /* current CSWB */
extern	double zc;	/* current SWB `borrow` */
extern	double zx;	/* SWB seed1 */
extern	double zy;	/* SWB seed2 */
extern	size_t qlen;/* length of Q array */


void NURBSurfacesMPS::InitiateRandNumGenerator (unsigned long x)
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
double NURBSurfacesMPS::RandNumGenerator()
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

void NURBSurfacesMPS::NURBSDartThrowing(size_t _num_active,size_t **_active_uw,size_t**_tmp_active_uw,double _s,size_t _nu,size_t _nw,size_t _n,//active cell/segemnts stuff
	                                  size_t&_num_samples,double**&_samples,double _r_input,double _tol,size_t _num_expected, // samples stuff
									  size_t* _post,size_t*_negt,// kd-tree stuff
	                                  double*** ctrl_p_ij,double _u_max,double _w_max,size_t *knot_u,size_t *knot_w,size_t K,size_t N,size_t L,size_t M,bool kts_u,bool kts_w)// nurb curve stuff
{
	//srand ( time(NULL) ); // activate for different experiments
	InitiateRandNumGenerator(rand());

	ConstructNURBS nurb;
	size_t num_darts,idart,rand_index,ii,jj,num_tmp_active,iactive;

	double RF(0.8),uu,ww,xx,yy,zz,uu_min,ww_min,uu_max,ww_max,ss;

	for(int iref=0;iref<30;iref++){
		cout<<"num_samples= "<<_num_samples<<endl;
		cout<<"num_active= "<<_num_active<<endl;
		num_darts=RF*_num_active;

		if(false){
			DrawActiveCells(_num_active,_active_uw,_s,_tol,
				ctrl_p_ij,_u_max,_w_max,knot_u,knot_w,K,N,L,M,kts_u,kts_w);
		}
			        

		for(idart=0;idart<num_darts;idart++){
			rand_index=size_t((_num_active-1)*RandNumGenerator());

			ii=_active_uw[rand_index][0];
			jj=_active_uw[rand_index][1];

			uu=(ii+RandNumGenerator())*_s;
			ww=(jj+RandNumGenerator())*_s;

			if(uu<0.0 || ww<0.0 || uu>_u_max || ww>_w_max){continue;}
			
			nurb.PointOnNURBSurface(uu,ww,N,M,K,L,ctrl_p_ij,knot_u,knot_w,xx,yy,zz,kts_u,kts_w,_tol,_u_max,_w_max);

			if(!Conflicting(xx,yy,zz,0,_samples,_post,_negt,_r_input)){
				//add sample if not conflicting
				_samples[_num_samples][0]=xx;
				_samples[_num_samples][1]=yy;
				_samples[_num_samples][2]=zz;
				_samples[_num_samples][3]=uu;
				_samples[_num_samples][4]=ww;
				_num_samples++;
				//add sample to kd-tree
				if(_num_samples>1){AddToKdTree(xx,yy,zz,_num_samples-1,_samples,_post,_negt);}

				if(double(_num_samples)-10>_num_expected){
					cout<<"Error (0) at NURBSDartThrowing().. _num_samples>_num_expected"<<endl;
					system("pause");
				}
			}
		}

		num_tmp_active=0;		

		for(iactive=0;iactive<_num_active;iactive++){
			ii=_active_uw[iactive][0];
			jj=_active_uw[iactive][1];

			uu_min=ii*_s;
			ww_min=jj*_s;
			uu_max=uu_min+_s;
			ww_max=ww_min+_s;
				  


			if(uu_min>=_u_max|| ww_max<0 || ww_min>=_w_max || uu_max<0){continue;}

			

			if(Covered(uu_min,ww_min,_s,0,_samples,_post,_negt,_r_input,_tol,
				      ctrl_p_ij,_u_max,_w_max,knot_u,knot_w,K,N,L,M,kts_u,kts_w)){continue;}

			_tmp_active_uw[num_tmp_active][0]=2*ii;
			_tmp_active_uw[num_tmp_active][1]=2*jj;
			num_tmp_active++;

			_tmp_active_uw[num_tmp_active][0]=2*ii+1;
			_tmp_active_uw[num_tmp_active][1]=2*jj;
			num_tmp_active++;

			_tmp_active_uw[num_tmp_active][0]=2*ii;
			_tmp_active_uw[num_tmp_active][1]=2*jj+1;
			num_tmp_active++;

			_tmp_active_uw[num_tmp_active][0]=2*ii+1;
			_tmp_active_uw[num_tmp_active][1]=2*jj+1;
			num_tmp_active++;

			if(double(num_tmp_active)-10 >_num_expected){
				cout<<"Error (1) at NURBSDartThrowing().. _num_active>_num_expected"<<endl;
				system("pause");
			}
		}
		if(num_tmp_active==0){break;}
		else{
			_num_active=num_tmp_active;
			for(iactive=0;iactive<num_tmp_active;iactive++){
				_active_uw[iactive][0]=_tmp_active_uw[iactive][0];
				_active_uw[iactive][1]=_tmp_active_uw[iactive][1];
			}
		}
		_s/=2.0;
	}
}

void NURBSurfacesMPS::AddToKdTree (double xx,double yy,double zz,size_t new_point,double**_samples,size_t*_post,size_t*_negt)
{
	// adding new point to the tree
	
	double dd;
	size_t k(0);

	while(true){
		
		if(k%3==0){dd=_samples[k][0]-xx;}
		else if(k>=1 && (k-1)%3==0){dd=_samples[k][1]-yy;}
		else if(k>=2 && (k-2)%3==0){dd=_samples[k][2]-zz;}
		else {
			cout<<"Error (0) at AddToKdTree()..!!"<<endl;
			system("pause");
		}
		
		if(dd<0){
			
			if(_post[k]==0){_post[k]=new_point;break;}
			else{k=_post[k]; continue;}
		}
		else{
			if(_negt[k]==0){_negt[k]=new_point;break;}
			else{k=_negt[k];continue;}
		}
	}
	
}
bool NURBSurfacesMPS::Conflicting(double xx, double yy, double zz,size_t k, double**_samples, size_t*_post,size_t*_negt,double _r_input)
{
	
	double dd;

	if(k%3==0){dd=_samples[k][0]-xx;}
	else if(k>=1 && (k-1)%3==0){dd=_samples[k][1]-yy;}
	else if(k>=2 && (k-2)%3==0){dd=_samples[k][2]-zz;}
	else {
		cout<<"Error (0) at Conflicting()..!!"<<endl;
		system("pause");
	}

	if(fabs(dd)<=3.0*_r_input){
		if(Dist(_samples[k][0],_samples[k][1],_samples[k][2],xx,yy,zz)<_r_input*_r_input){
			return true;
		}

		if(_post[k]==0 && _negt[k]==0){return false;}
		if(_post[k]!=0 && _negt[k]==0){return Conflicting(xx,yy,zz,_post[k],_samples,_post,_negt,_r_input);}
		if(_post[k]==0 && _negt[k]!=0){return Conflicting(xx,yy,zz,_negt[k],_samples,_post,_negt,_r_input);}
		if(_post[k]!=0 && _negt[k]!=0){

			if(!Conflicting(xx,yy,zz,_negt[k],_samples,_post,_negt,_r_input)){
				if(!Conflicting(xx,yy,zz,_post[k],_samples,_post,_negt,_r_input)){return false;}
				else{return true;}
			}
			else {return true;}
		}
	}
	else if(dd>=0.0 && _negt[k]!=0) {return Conflicting(xx,yy,zz,_negt[k],_samples,_post,_negt,_r_input);}
	else if(dd< 0.0 && _post[k]!=0) {return Conflicting(xx,yy,zz,_post[k],_samples,_post,_negt,_r_input);}
	else {return false;}


}
double NURBSurfacesMPS::Dist(double x1,double y1, double z1, double x2, double y2, double z2)
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
bool NURBSurfacesMPS::Covered(double uu_min,double ww_min,double _s,size_t k,double**_samples,size_t*_post,size_t*_negt,double _r_input,double _tol,
	                          double*** ctrl_p_ij,double _u_max,double _w_max,size_t *knot_u,size_t *knot_w,size_t K,size_t N,size_t L,size_t M,bool kts_u,bool kts_w)
{
	ConstructNURBS nurb;
	double xx1,yy1,zz1,dd;
	nurb.PointOnNURBSurface(uu_min,ww_min,N,M,K,L,ctrl_p_ij,knot_u,knot_w,xx1,yy1,zz1,kts_u,kts_w,_tol,_u_max,_w_max);

	if(k%3==0){dd=_samples[k][0]-xx1;}
	else if(k>=1 && (k-1)%3==0){dd=_samples[k][1]-yy1;}
	else if(k>=2 && (k-2)%3==0){dd=_samples[k][2]-zz1;}
	else {
		cout<<"Error (0) at Covered()..!!"<<endl;
		system("pause");
	}

	if(fabs(dd)<=3.0*_r_input){

		if(Dist(_samples[k][0],_samples[k][1],_samples[k][2],xx1,yy1,zz1)<_r_input*_r_input){
			double xx2,yy2,zz2;
			nurb.PointOnNURBSurface(uu_min+_s,ww_min,N,M,K,L,ctrl_p_ij,knot_u,knot_w,xx2,yy2,zz2,kts_u,kts_w,_tol,_u_max,_w_max);
			if(Dist(_samples[k][0],_samples[k][1],_samples[k][2],xx2,yy2,zz2)<_r_input*_r_input){
				double xx3,yy3,zz3;
				nurb.PointOnNURBSurface(uu_min,ww_min+_s,N,M,K,L,ctrl_p_ij,knot_u,knot_w,xx3,yy3,zz3,kts_u,kts_w,_tol,_u_max,_w_max);
				if(Dist(_samples[k][0],_samples[k][1],_samples[k][2],xx3,yy3,zz3)<_r_input*_r_input){
					double xx4,yy4,zz4;
					nurb.PointOnNURBSurface(uu_min+_s,ww_min+_s,N,M,K,L,ctrl_p_ij,knot_u,knot_w,xx4,yy4,zz4,kts_u,kts_w,_tol,_u_max,_w_max);
					if(Dist(_samples[k][0],_samples[k][1],_samples[k][2],xx4,yy4,zz4)<_r_input*_r_input){
						return true;
					}
				}
			}
		}

		if(_post[k]==0 && _negt[k]==0){return false;}
		if(_post[k]!=0 && _negt[k]==0){return Covered(uu_min,ww_min,_s,_post[k],_samples,_post,_negt,_r_input,_tol, 
			                                          ctrl_p_ij,_u_max,_w_max,knot_u,knot_w,K,N,L,M,kts_u,kts_w);}
		if(_post[k]==0 && _negt[k]!=0){return Covered(uu_min,ww_min,_s,_negt[k],_samples,_post,_negt,_r_input,_tol,
			                                          ctrl_p_ij,_u_max,_w_max,knot_u,knot_w,K,N,L,M,kts_u,kts_w);}
		if(_post[k]!=0 && _negt[k]!=0){
			if(Covered(uu_min,ww_min,_s,_negt[k],_samples,_post,_negt,_r_input,_tol, 
				       ctrl_p_ij,_u_max,_w_max,knot_u,knot_w,K,N,L,M,kts_u,kts_w)){return true;}
			else if(Covered(uu_min,ww_min,_s,_post[k],_samples,_post,_negt,_r_input,_tol, 
				       ctrl_p_ij,_u_max,_w_max,knot_u,knot_w,K,N,L,M,kts_u,kts_w)){return true;}
			else {return false;}				
		}
	}

	else if(dd>=0.0 && _negt[k]!=0) {return Covered(uu_min,ww_min,_s,_negt[k],_samples,_post,_negt,_r_input,_tol, 
		                                            ctrl_p_ij,_u_max,_w_max,knot_u,knot_w,K,N,L,M,kts_u,kts_w);}
	else if(dd< 0.0 && _post[k]!=0) {return Covered(uu_min,ww_min,_s,_post[k],_samples,_post,_negt,_r_input,_tol,  
		                                            ctrl_p_ij,_u_max,_w_max,knot_u,knot_w,K,N,L,M,kts_u,kts_w);}
	else {return false;}
}
void NURBSurfacesMPS::DrawActiveCells(size_t num_active, size_t** _active_uw, double _s,double _tol,
	                 double*** ctrl_p_ij,double _u_max,double _w_max,size_t *knot_u,size_t *knot_w,size_t K,size_t N,size_t L,size_t M,bool kts_u,bool kts_w)// nurb curve stuff)
{
	ConstructNURBS nurb;
	size_t iactive,ii,jj,num(0);
	double uu_min,ww_min,uu_max,ww_max,x,y,z;
	fstream file("active_cells_surface.obj", ios::out);
	file.precision(30);


	for(iactive=0;iactive<num_active;iactive++){
			ii=_active_uw[iactive][0];
			jj=_active_uw[iactive][1];

			uu_min=ii*_s;
			ww_min=jj*_s;
			uu_max=uu_min+_s;
			ww_max=ww_min+_s;

			
			if(iactive==49){
				size_t dsgf=45;
			}

			nurb.PointOnNURBSurface(uu_min,ww_min,N,M,K,L,ctrl_p_ij,knot_u,knot_w,x,y,z,kts_u,kts_w,_tol,_u_max,_w_max);
			file<<"v "<<x<<" "<<y<<" "<<z<<endl;

			if(! (x<100 && x>-100 &&  y<100 && y>-100 && z<100 && z>-100)){
				size_t sfgs=45;
			}

			nurb.PointOnNURBSurface(uu_min,ww_max,N,M,K,L,ctrl_p_ij,knot_u,knot_w,x,y,z,kts_u,kts_w,_tol,_u_max,_w_max);
			file<<"v "<<x<<" "<<y<<" "<<z<<endl;

			nurb.PointOnNURBSurface(uu_max,ww_max,N,M,K,L,ctrl_p_ij,knot_u,knot_w,x,y,z,kts_u,kts_w,_tol,_u_max,_w_max);
			file<<"v "<<x<<" "<<y<<" "<<z<<endl;

			nurb.PointOnNURBSurface(uu_max,ww_min,N,M,K,L,ctrl_p_ij,knot_u,knot_w,x,y,z,kts_u,kts_w,_tol,_u_max,_w_max);
			file<<"v "<<x<<" "<<y<<" "<<z<<endl;



			num++;
	}

	size_t V,i(1);
	for(V=0;V<num;V++){
		file<<"f "<<i<<" "<<i+1<<" "<<i+2<<" "<<i+3<<endl;		
		i+=4;
	}

}