class ProjectPointToNURBCurve
{
public:

	double MinDist(double,double,double&,double**,double,size_t*,size_t,size_t,bool);
private:
	double Dist(double,double,double,double,double,double);
	void PlotThatPoint(double,double,double,double,double**,size_t);
};