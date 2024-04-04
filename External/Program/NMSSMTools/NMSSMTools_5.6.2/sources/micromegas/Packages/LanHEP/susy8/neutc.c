

double neut1_(double *MZ,double *MW, double *SW, double *MG1, 
	  	double *MG2, double *mue, double *tb);
		
double neut2_(double *tk, double *o1, double *o2, double *o3);

double neut1(double MZ,double MW, double SW, double MG1, 
	  	double MG2, double mue, double tb)
	{
		return neut1_(&MZ, &MW, &SW, &MG1, &MG2, &mue, &tb);
	}
	
double neut2(double tk, double o1, double o2, double o3)
	{
		return neut2_(&tk, &o1, &o2, &o3);
	}
