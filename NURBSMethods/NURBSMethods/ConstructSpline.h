class ConstructSpline
{

public:	
	void PointOnSplineCurve(double,size_t,size_t,double**,size_t*&,double&,double&,bool,double,double);	
	void PlotSpline_ps(size_t,size_t,double**,size_t*,bool,double,double);
	void PointOnSplineSurface(double,double,size_t,size_t,size_t,size_t,double***,size_t*&,size_t*&,double&,double&,double&,bool,bool,double,double,double);
	void PlotSpline_surface(size_t,size_t,size_t,size_t,double***,size_t*,size_t*,bool,bool,double,double,double);
private:
	void KnotVector(size_t, size_t,size_t*&);
	size_t N_i_1(double,size_t,size_t*);
	double BasisFunction(size_t,size_t,size_t,double,size_t*);
};