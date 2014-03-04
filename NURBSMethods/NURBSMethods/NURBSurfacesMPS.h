class NURBSurfacesMPS
{

public:
	void NURBSDartThrowing(size_t,size_t**,size_t**,double,size_t,size_t,size_t,//active cell/segemnts stuff
	                                  size_t&,double**&,double,double,size_t, // samples stuff
									  size_t*,size_t*,// kd-tree stuff
	                                  double***,double,double,size_t*,size_t*,size_t,size_t,size_t,size_t,bool,bool);// nurb curve stuff

private:
	double RandNumGenerator();
	void InitiateRandNumGenerator(unsigned long);
	void AddToKdTree(double,double,double,size_t,double**,size_t*,size_t*);
	bool Conflicting(double,double,double,size_t,double**,size_t*,size_t*,double);
	double Dist(double,double,double,double,double,double);
	bool Covered(double,double,double,size_t,double**,size_t*,size_t*,double,double,
            	 double***,double,double,size_t*,size_t*,size_t,size_t,size_t,size_t,bool,bool);	
	void DrawActiveCells(size_t,size_t**,double,double,
		                 double***,double,double,size_t*,size_t*,size_t,size_t,size_t,size_t,bool,bool);// nurb curve stuff);


};