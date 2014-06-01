#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <string.h>  
#include <fstream>
#include <sstream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <bitset>
#include <time.h>
#include <direct.h>
using namespace std;


double **_surface;// surface points coordinates
size_t _num_nodes;// number of points of the input surafce mesh
size_t _num_tri; // number of triangle of the input surface mesh 
size_t ** _surface_tri;// connectivity of the input surface
                       //legnth= _num_nodes
                       //each row:
                       //#n q0 q1 q2 q3 qn-1
                       //n=number of points connected to 
size_t**_surface_tri2;// connectivity of the input surface
                      //length=_num_tri
                      //each row 
                      //t1 t2 t3 n1 n2 n3
                      //t's are the triangle three vertics, n's are the triangles three neighbour triangles
                      //for "bad" input mesh, deplicated neighbours maybe found 
size_t _num_layers;//number of layers of triangles to propoagte through 

size_t *_neighbor_tri,*_neighbor_tri2; //stores the neighbour triangles in number of layers = _num_layers 

inline double Dot(double xv1, double yv1, double zv1, double xv2, double yv2, double zv2)
{
	double dot;
	dot=xv1*xv2 + yv1*yv2 + zv1*zv2;

	return dot;
}
inline double PointTriangleDistance(size_t tri, double xp, double yp, double zp,double&x_new,double&y_new,double&z_new)
{
   // http://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
   // http://www.mathworks.com/matlabcentral/fileexchange/22857-distance-between-a-point-and-a-triangle-in-3d/content/pointTriangleDistance.m
	//        ^t
	//  \     |
	//   \reg2|
	//    \   |
	//     \  |
	//      \ |
	//       \|
	//        *P2
	//        |\
	//        | \
	//  reg3  |  \ reg1
	//        |   \
	//        |reg0\ 
	//        |     \ 
	//        |      \ P1
	// -------*-------*------->s
	//        |P0      \ 
	//  reg4  | reg5    \ reg6

	size_t p1(_surface_tri2[tri][0]),p2(_surface_tri2[tri][1]),p3(_surface_tri2[tri][2]);
	double xp1(_surface[p1][0]),yp1(_surface[p1][1]),zp1(_surface[p1][2]),
      	   xp2(_surface[p2][0]),yp2(_surface[p2][1]),zp2(_surface[p2][2]),
		   xp3(_surface[p3][0]),yp3(_surface[p3][1]),zp3(_surface[p3][2]);

	double xE0(xp2-xp1),yE0(yp2-yp1),zE0(zp2-zp1),
           xE1(xp3-xp1),yE1(yp3-yp1),zE1(zp3-zp1),
		   xD(xp1-xp),yD(yp1-yp),zD(zp1-zp);

	double a=Dot(xE0,yE0,zE0,xE0,yE0,zE0);
	double b=Dot(xE0,yE0,zE0,xE1,yE1,zE1);
	double c=Dot(xE1,yE1,zE1,xE1,yE1,zE1);
	double d=Dot(xE0,yE0,zE0,xD,yD,zD);
	double e=Dot(xE1,yE1,zE1,xD,yD,zD);
	double f=Dot(xD,yD,zD,xD,yD,zD);

	double det=a*c-b*b;
	double s=b*e-c*d;
	double t=b*d-a*e;
	double invDet,dist,numer,denom,tmp0,tmp1;
	//Terible tree of conditionals to determine in which region of the diagram
	// shown above the projection of the point into the triangle-plane lies.

  if(s+t<=det){
	  if(s<0){
		  if(t<0){
			  //region 4
			   if (d < 0){
				   t=0;
				    if (-1.0*d >= a){
						s=1;
						dist=a + 2.0*d + f;
					}
					else{
						s = -1.0*d/a;
						dist = d*s + f;
					}
			   }
			   else{
				   s = 0;
				   if (e >= 0){
					   t = 0;
					   dist = f;
				   }
				   else{
					   if (-e >= c){
						   t = 1;
						   dist = c + 2.0*e + f;
					   }
					   else{
						   t = -1.0*e/c;
						   dist = e*t + f;
					   }

				   }

			   }
			  
		  }
		  else
		  {
			  //region 3
			  s=0;
			  if (e >= 0){
				  t = 0;
				  dist = f;
			  }
			  else{
				  if (-1.0*e >= c){
					  t=1;
					  dist=c+2.0*e+f;
				  }
				  else{
					  t = -1.0*e/c;
					  dist = e*t + f;
				  }
			  }

		  }
	  }
	  else if(t<0){
		  //region 5
		  t=0;
		  if (d >= 0){
			  s = 0;
			  dist = f;
		  }
		  else{
			  if (-1.0*d >= a){
				  s=1;
				  dist=a+2.0*d+f;
			  }
			  else{
				  s = -1.0*d/a;
				  dist = d*s + f;
			  }
		  }

	  }
	  else{
		  //region 0
		  invDet=1.0/det;
		  s=s*invDet;
		  t=t*invDet;
		  dist=s*(a*s + b*t + 2.0*d) +
	          t*(b*s + c*t + 2.0*e) + f;
	  }
  }
  else{
	  if(s<0){
		  //region 2
		  tmp0 = b + d;
		  tmp1 = c + e;
		   if (tmp1 > tmp0){ // minimum on edge s+t=1
			   numer = tmp1 - tmp0;
			   denom = a - 2.0*b + c;
			   if(numer >= denom){
				   s = 1;
				   t = 0;
				   dist = a + 2.0*d + f;
			   }
			   else{
				   s = numer/denom;
				   t = 1.0-s;
				   dist = s*(a*s + b*t + 2.0*d)
					            + t*(b*s + c*t + 2*e) + f;
			   }
		   }
		   else{
			    s = 0;
				if(tmp1 <= 0){
					t = 1;
					dist= c + 2.0*e + f;
				}
				else{
					if (e >= 0){
						t = 0;
						dist = f;
					}
					else{
						t = -1.0*e/c;
						dist = e*t + f;
					}
				}
		   }
	  }
	  else if(t<0){
		  //region 6

		  tmp0 = b + e;
		  tmp1 = a + d;
		   if (tmp1 > tmp0){ // minimum on edge s+t=1
			   numer = tmp1 - tmp0;
			   denom = a - 2.0*b + c;
			   if(numer >= denom){
				   t = 1;
				   s = 0;
				   dist = c + 2.0*e + f;
			   }
			   else{
				   t = numer/denom;
				   s = 1.0-t;
				   dist = s*(a*s + b*t + 2.0*d)
                      + t*(b*s + c*t + 2.0*e) + f;
			   }
		   }
		   else{
			    t = 0;
				if(tmp1 <= 0){
					s = 1;
					dist= a + 2.0*d + f;
				}
				else{
					if (d >= 0){
						s = 0;
						dist = f;
					}
					else{
						s = -d/a;
						dist = d*s + f;
					}
				}
		   }
	  }
	  else{
		  //region 1
		  numer=c+e-b-d;
		  if(numer<=0){
			  s=0;
			  t=1;
			  dist=c + 2.0*e+f;
		  }
		  else{
			  denom=a - 2.0*b + c;
			   if (numer >= denom){
				   s = 1;
				   t = 0;
				   dist= a + 2.0*d + f;
			   }
			   else{
				   s = numer/denom;
				   t = 1-s;
				   dist = s*(a*s + b*t + 2.0*d)+
					             t*(b*s + c*t + 2.0*e) + f;
			   }
		  }
	  }
  }

  x_new=xp1+s*xE0+t*xE1;
  y_new=yp1+s*yE0+t*yE1;
  z_new=zp1+s*zE0+t*zE1;
   return sqrt(dist);

}
void GetTriNeighbors(size_t tri_num)
{
	size_t t1,t2,t3,V,tt;
	t1=_surface_tri2[tri_num][3];
	t2=_surface_tri2[tri_num][4];
	t3=_surface_tri2[tri_num][5];
	bool in=true;

	for(V=1;V<=_neighbor_tri[0];V++){
		tt=_neighbor_tri[V];
		
		if(t1==tt){in=false;break;}
	}
	if(in){
		_neighbor_tri[0]++;
		_neighbor_tri[_neighbor_tri[0]]=t1;
	}
	in=true;


	
	for(V=1;V<=_neighbor_tri[0];V++){
		tt=_neighbor_tri[V];
		
		if(t2==tt){in=false;break;}
	}
	if(in){
		_neighbor_tri[0]++;
		_neighbor_tri[_neighbor_tri[0]]=t2;
	}
	in=true;
	

	
	for(V=1;V<=_neighbor_tri[0];V++){
		tt=_neighbor_tri[V];
		
		if(t3==tt){in=false;break;}
	}
	if(in){
		_neighbor_tri[0]++;
		_neighbor_tri[_neighbor_tri[0]]=t3;
	}
	in=true;

}
inline void ProjectPointToSurface(double &xx, double&yy, double&zz, size_t&t1)
{
	//projecting point (xx,yy,zz) to the surface 
	//where t1 is a rough estimation to the closest triangle before projection
	//return xx,yy,zz as the projection and t1 as the triangle containing the point after projection
	//if t1 is unknown (in very rare cases) we shall loop over all triangles in the surface (requires increasing _num_layers maybe)
	size_t V,d,counter,tri_replace,start,num;
	double dist(10E6),dd,x_replace,y_replace,z_replace,x2,y2,z2;
		
	_neighbor_tri[0]=4;	/*collecting neighbour triangles*/
	_neighbor_tri[1]=t1;
	_neighbor_tri[2]=_surface_tri2[t1][3];
	_neighbor_tri[3]=_surface_tri2[t1][4];
	_neighbor_tri[4]=_surface_tri2[t1][5];
	
	counter=0;
	
	start=1;
	while(counter<_num_layers){
		num=_neighbor_tri[0];
		for(d=start;d<=num;d++){
			
			GetTriNeighbors(_neighbor_tri[d]);
		}
		start=num;
		counter++;
	}
	

	for(V=1;V<=_neighbor_tri[0];V++){/*calculate min distance*/

		dd=PointTriangleDistance(_neighbor_tri[V],xx,yy,zz,x2,y2,z2);
		
		if(dd<dist){
			dist=dd;
			x_replace=x2;
			y_replace=y2;
			z_replace=z2;
			tri_replace=_neighbor_tri[V];
		}		
	}



	xx=x_replace; /*new location*/
	yy=y_replace;
	zz=z_replace;


	t1=tri_replace;

		
	
}

inline void SortSurface(size_t ip)
{
	//for surface read
	size_t point_1, point_2,tmp,next_point,jpoint;
	double tmpv;
	bool Connectivity=false;
	

	for ( jpoint=1; jpoint<_surface_tri[ip][0];jpoint++){

		point_1=_surface_tri[ip][jpoint]; // leaving the first entry poitn as it is

		for ( next_point=jpoint+1; next_point<=_surface_tri[ip][0];next_point++){// looping arround the rest of the nodes and check the connectivity 
		
			point_2=_surface_tri[ip][next_point];
			Connectivity=false;

			 //check on the connectivity between jpoint and next_point
			
			for (size_t R=1; R<=_surface_tri[point_1][0]; R++){
				if (_surface_tri[point_1][R]==point_2){Connectivity=true;break;}				
			}
			for (size_t R=1; R<=_surface_tri[point_2][0]; R++){
				if (_surface_tri[point_2][R]==point_1){Connectivity=true;break;}				
			}


			if (Connectivity==true){  // swapping if connected 
			    
				tmp=_surface_tri[ip][jpoint+1];
				_surface_tri[ip][jpoint+1]=point_2;
				_surface_tri[ip][next_point]=tmp;
				

				break;
			}
		}
	}
	
}
inline bool IsThere(size_t inode, size_t ipoint)
{
	//for surface read
	for(size_t V=1; V<=_surface_tri[inode][0];V++){

		if(_surface_tri[inode][V]==ipoint){return true;}
	}

	return false;
}
inline void Plot()
{
	size_t V,d;

	fstream file55("surface.obj",ios::out);		
	for(V=0; V<_num_nodes;V++){			
		file55<<"v "<<_surface[V][0]<<" "<<_surface[V][1]<<" "<<_surface[V][2]<<endl;
	}
	
	for(V=0;V<_num_nodes;V++){
		for(d=1; d<=_surface_tri[V][0];d++){
			
			if(d==_surface_tri[V][0]){
				file55<<"f "<<V+1<<" "<<_surface_tri[V][d]+1<<" "<<_surface_tri[V][1]+1<<endl;
			}else{
				file55<<"f "<<V+1<<" "<<_surface_tri[V][d]+1<<" "<<_surface_tri[V][d+1]+1<<endl;
			}
		}
	}
	
}
int main(int argc, char **argv)
{
	ifstream input;
	input.open("Sphere2.off");	

	#pragma region ReadSurface 
	//*******Reading and setting initial values starts here******//
	size_t V;
	input>>_num_nodes;
	input>>_num_tri;
	
	_surface = new double* [_num_nodes];
	_surface_tri=new size_t*[_num_nodes];
	_surface_tri2=new size_t*[_num_tri];
	size_t** temp_node_ele=new size_t*[_num_nodes]; // temp array to store triangles belong to each node to find out each element neighbors

	for(V=0;V<_num_nodes;V++){
		_surface[V]=new double [3];
		input>>_surface[V][0]>>_surface[V][1]>>_surface[V][2]; //read samples coordinates

		temp_node_ele[V]=new size_t[50];
		temp_node_ele[V][0]=0;	

		_surface_tri[V]=new size_t[30+1];
		_surface_tri[V][0]=0;
	}
		
	size_t node1,node2,node3,tuna;
	double dx,dy,dz,dist_sq;

	double max(0),min(500000000000);

	
	for (V=0;V<_num_tri;V++){ // read delaunay triangles 
		input>>tuna>>node1>>node2>>node3;	
		node1-=1;
		node2-=1;
		node3-=1;
		
		_surface_tri2[V]=new size_t[6];

		_surface_tri2[V][0]=node1;
		_surface_tri2[V][1]=node2;
		_surface_tri2[V][2]=node3;

		temp_node_ele[node1][0]+=1;
		temp_node_ele[node1][temp_node_ele[node1][0]]=V;

		temp_node_ele[node2][0]+=1;
		temp_node_ele[node2][temp_node_ele[node2][0]]=V;

		temp_node_ele[node3][0]+=1;
		temp_node_ele[node3][temp_node_ele[node3][0]]=V;
		
			
		
		if(!IsThere(node1,node2)){_surface_tri[node1][0]++; _surface_tri[node1][_surface_tri[node1][0]]=node2;}
		if(!IsThere(node1,node3)){_surface_tri[node1][0]++; _surface_tri[node1][_surface_tri[node1][0]]=node3;}

		if(!IsThere(node2,node1)){_surface_tri[node2][0]++; _surface_tri[node2][_surface_tri[node2][0]]=node1;}
		if(!IsThere(node2,node3)){_surface_tri[node2][0]++; _surface_tri[node2][_surface_tri[node2][0]]=node3;}

		if(!IsThere(node3,node1)){_surface_tri[node3][0]++; _surface_tri[node3][_surface_tri[node3][0]]=node1;}
		if(!IsThere(node3,node2)){_surface_tri[node3][0]++; _surface_tri[node3][_surface_tri[node3][0]]=node2;}
	}
		
	
	for(V=0;V<_num_nodes;V++){
		SortSurface(V);
	}
	

	size_t p1,p2,p3,d;
	int p;
	bool finish=false;
	size_t t1,t2;

	for (V=0; V<_num_tri; V++) // finiding the neighbors by looping over each node's triangles
	{
		p1=_surface_tri2[V][0]; p2=_surface_tri2[V][1]; p3=_surface_tri2[V][2];

		for(d=1; d<=temp_node_ele[p1][0]; d++)
		{
			t1=temp_node_ele[p1][d];
			if(t1==V){continue;}

			for(p=1;p<=temp_node_ele[p2][0];p++)
			{
				t2=temp_node_ele[p2][p];
				if(t2==V){continue;}

				if(t1==t2)
				{
					_surface_tri2[V][3]=t1;
					finish=true;
					break;
				}
			}

			if(finish){break;}
		}

		
		finish=false;



		for(d=1; d<=temp_node_ele[p2][0]; d++)
		{
			t1=temp_node_ele[p2][d];
			if(t1==V){continue;}

			for(p=1;p<=temp_node_ele[p3][0];p++)
			{
				t2=temp_node_ele[p3][p];
				if(t2==V){continue;}

				if(t1==t2)
				{
					_surface_tri2[V][4]=t1;
					finish=true;
					break;
				}
			}

			if(finish){break;}
		}

		finish=false;



		for(d=1; d<=temp_node_ele[p3][0]; d++)
		{
			t1=temp_node_ele[p3][d];
			if(t1==V){continue;}

			for(p=1;p<=temp_node_ele[p1][0];p++)
			{
				t2=temp_node_ele[p1][p];
				if(t2==V){continue;}

				if(t1==t2)
				{
					_surface_tri2[V][5]=t1;
					finish=true;
					break;
				}
			}

			if(finish){break;}
		}

		finish=false;

	}	
#pragma endregion
	



	Plot(); //input surface

	_neighbor_tri=new size_t[5000+1];
	_neighbor_tri[0]=0;


	//generate a point 
	//srand ( time(NULL) ); /*Random seed*/
	double x,y,z;
	size_t t=0;
	x=(double)rand() / RAND_MAX;
	y=(double)rand() / RAND_MAX;
	z=(double)rand() / RAND_MAX;

	x+=0.5; //scale
	y+=0.5;
	z+=0.5;

	_num_layers=50; //should change from model to another 

	cout<<"Point before projection ("<<x<<","<<y<<","<<z<<")"<<endl;
	//we pass (x,y,z,t) of the point to be projected
	//where (t) is a rough estimation to the closest triangle to the point to be projected
	// we collect the neighbour triangles around (t) in number of layers=_num_layers
	// then we calculate the min distance between the point and each of these triangles 
	// we return the closest point and the triangle contain this point 

	//if (t) is unknown (in very rare cases) we shall loop over all triangles in the surface (requires increasing _num_layers maybe)
	ProjectPointToSurface(x,y,z,t);
	cout<<"\nPoint after projection ("<<x<<","<<y<<","<<z<<")"<<endl;

	
	
	cout<<"\n Press Any Key To Exit..."<<endl;
	cin.get();
		
	return 0;

}