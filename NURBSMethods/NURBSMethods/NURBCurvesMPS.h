class NURBCurvesMPS
{
public:
	void NURBSDartThrowing(size_t,double*,double,double*,  //active cell/segemnts stuff
	                       size_t&,double**&,double,double, // samples stuff
	                       double**,double,size_t*,size_t,size_t,bool);// nurb curve stuff


private:
	void InitiateRandNumGenerator (unsigned long);
	double RandNumGenerator();
	bool Covered(double,double,double,double,double,double,double**,size_t,double,double**,double,size_t*,size_t,size_t,bool,double);
	double DistOnNurb(double,double,double**,double,size_t*,size_t,size_t,bool,double);
	double Dist(double,double,double,double,double,double);
	void PlotThatPoint(double,double,size_t,double**,size_t,double);
	bool Conflicting(double,double,double,size_t,double**,double,double**,double,size_t*,size_t,size_t,bool,double);


};