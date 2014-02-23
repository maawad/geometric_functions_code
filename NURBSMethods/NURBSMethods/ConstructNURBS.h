class ConstructNURBS
{

public:	
	void PointOnNURBCurve(double,size_t,size_t,double**,size_t*&,double&,double&,bool);	
	void PlotNURB_ps(size_t,size_t,double**,size_t*,bool,double);
	void PointOnNURBSurface(double,double,size_t,size_t,size_t,size_t,double***,size_t*&,size_t*&,double&,double&,double&,bool,bool);
	void PlotNURB_surface(size_t,size_t,size_t,size_t,double***,size_t*,size_t*,bool,bool,double,double);
private:
	void KnotVector(size_t, size_t,size_t*&);
	size_t N_i_1(double,size_t,size_t*);
	double BasisFunction(size_t,size_t,size_t,double,size_t*);
};