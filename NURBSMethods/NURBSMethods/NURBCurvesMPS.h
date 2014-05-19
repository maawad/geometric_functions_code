class NURBCurvesMPS
{
public:
	void NURBSDartThrowing(size_t,double*,double,double*,  //active cell/segemnts stuff
	                       size_t&,double**&,double,double, // samples stuff
	                       double**,double,size_t*,size_t,size_t,bool,// nurb curve stuff
						   size_t);


private:
	void InitiateRandNumGenerator (unsigned long);
	double RandNumGenerator();
	bool Covered(double,double,double,double,double,double,double**,size_t,double,double**,double,size_t*,size_t,size_t,bool,double,size_t);
	double DistOnNurb(double,double,double**,double,size_t*,size_t,size_t,bool,double);
	double Dist(double,double,double,double,double,double);
	void PlotThatPoint(double,double,size_t,double**,size_t,double);
	bool Conflicting(double,double,double,size_t,double**,double,double**,double,size_t*,size_t,size_t,bool,double,size_t);
	double DistOnSpline(double u1, double u2,double** ctrl_p,double _u_max,size_t *knot,size_t K,size_t N,bool kts,double _tol);


};