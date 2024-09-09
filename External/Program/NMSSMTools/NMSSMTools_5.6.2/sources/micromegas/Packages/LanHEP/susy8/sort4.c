#include <math.h>
#include <stdio.h>


double 
sort4 (double m1, double m2, double m3, double m4, double dn)
{
  int n;
  int i, f = 0;
  double m[4];
  n = (int) floor (dn - 0.5);
  m[0] = m1;
  m[1] = m2;
  m[2] = m3;
  m[3] = m4;

  do
    {
      f = 0;
      for (i = 0; i < 3; i++)
	if (fabs (m[i]) > fabs (m[i + 1]))
	  {
	    double tmp;
	    tmp = m[i];
	    m[i] = m[i + 1];
	    m[i + 1] = tmp;
	    f = 1;
	  }

    }
  while (f);

  return m[n];
}
