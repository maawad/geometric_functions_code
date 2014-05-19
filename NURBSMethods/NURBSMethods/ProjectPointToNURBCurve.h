class ProjectPointToNURBCurve
{
public:

	double MinDist(double,double,double&,double**,double,size_t*,size_t,size_t,bool,double,size_t);
private:
	double Dist(double,double,double,double,double,double);
	void PlotThatPoint(double,double,double,double,double**,size_t);
};