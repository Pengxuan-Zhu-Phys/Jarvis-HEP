
#include <math.h>

double neutcp1_(double *MZ,double *MW, double *SW, double *MG1, 
	  	double *MG2, double *mue, double *tb, double *th1, double *th2, double
		*thmu);
		
double neutcp2_(double *tk, double *o1, double *o2, double *o3);

double neutcp1(double MZ,double MW, double SW, double MG1, 
	  	double MG2, double mue, double tb, double th1, double th2, double thmu)
	{
		return neutcp1_(&MZ, &MW, &SW, &MG1, &MG2, &mue, &tb, &th1, &th2, &thmu);
	}
	
double neutcp2(double tk, double o1, double o2, double o3)
	{
		return neutcp2_(&tk, &o1, &o2, &o3);
	}

double chacp1_(double *MZ,double *MW, double *SW,  
	  	double *MG2, double *mue, double *tb, double *th2, double
		*thmu);
		
double chacp2_(double *tk, double *o1, double *o2, double *o3);


double chacp1(double MZ,double MW, double SW, 
	  	double MG2, double mue, double tb, double th2, double thmu)
	{
		return chacp1_(&MZ, &MW, &SW, &MG2, &mue, &tb, &th2, &thmu);
	}
	
double chacp2(double tk, double o1, double o2, double o3)
	{
		return chacp2_(&tk, &o1, &o2, &o3);
	}
