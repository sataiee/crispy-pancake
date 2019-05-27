#include <math.h>
#include "mp.h"

#define power1 4.44444444e-2
#define power2 2.381e-2
#define power3 2.267e-1

#define t234 1.6e3
#define t456 5.7e3
#define t678 2.28e6

/* coefficients for opacity laws 1, 2, and 3 in cgs units */
#define ak1 2.e-4
#define ak2 2.e16
#define ak3 5.e-3

/* coefficients for opacity laws 3, 4, 5, 6, 7, and 8 in T_4 units */
#define bk3 50.0
#define bk4 2.e-2
#define bk5 2.e4
#define bk6 1.e4
#define bk7 1.5e10
#define bk8 0.348

double lin(double temperature, double density)
{
  /* test T against (T_23 * T_34 * T_34)^(1/3) */
  if (temperature > t234*pow(density, power1)) 
    {
    /* to avoid overflow */
    double ts4 = 1.0e-4*temperature;
    double density13 = pow(density, 1.0/3.0);
    double density23 = density13*density13;
    double ts42 = ts4 * ts4;
    double ts44 = ts42 * ts42;
    double ts48 = ts44 * ts44;

    /* test T against (T_45 * T_56)^0.5 */
    if (temperature > t456*pow(density, power2)) 
        {
      if ((temperature < t678*pow(density, power3)) || 
                (density <= 1.0e-10)) 
            {
        /* disjoint opacity laws for 5, 6, and 7 */
        double o5 = bk5*density23*ts42*ts4;
        double o6 = bk6*density13*ts48*ts42;
        double o7 = bk7*density/(ts42*sqrt(ts4));

        /* parameters used for smoothing */
        double o6an = o6*o6;
        double o7an = o7*o7;

        /* smoothed and continuous opacity law for regions 5, 6, and 7 */
        return pow( pow(o6an*o7an/(o6an+o7an), 2)
                           +pow(o5/(1.0+pow(ts4/(1.1*pow(density,0.04762)),10.0)),4.0),0.25);
      } 
            else 
            {
        /* disjoint opacity laws for 7 and 8 */
        double o7 = bk7*density/(ts42*sqrt(ts4));
        double o8 = bk8;

        /* parameters used for smoothing */
        double o7an = o7*o7;
        double o8an = o8*o8;

        /* smoothed and continuous opacity law for regions 7 and 8 */
        return pow(o7an*o7an+o8an*o8an,0.25);

        /* no scattering */
        //return bk7*density/(ts42*sqrt(ts4));
      }
    } 
        else 
        {
      /*  disjoint opacity laws for 3, 4, and 5 */
      double o3 = bk3*ts4;
      double o4 = bk4*density23/(ts48*ts4);
      double o5 = bk5*density23*ts42*ts4;

      /* parameters used for smoothing */
      double o4an = o4*o4*o4*o4;
      double o3an = o3*o3*o3*o3;

      /* smoothed and continuous opacity law for regions 3, 4, and 5 */
      return pow( (o4an*o3an/(o4an+o3an))
                       +pow(o5/(1.0+6.561e-5/ts48),4),0.25);
    }
  } 
    else 
    {
    /* different powers of temperature */
    double t2 = temperature*temperature;
    double t4 = t2*t2;
    double t8 = t4*t4;
    double t10 = t8*t2;

    /* disjoint opacity laws */
    double o1 = ak1*t2;
    double o2 = ak2*temperature/t8;
    double o3 = ak3*temperature;

    /* parameters used for smoothing */
    double o1an = o1*o1;
    double o2an = o2*o2;

    /* smoothed and continuous opacity law for regions 1, 2, and 3 */
    return pow(pow(o1an*o2an/(o1an+o2an),2)+pow(o3/(1+1.e22/t10),4),0.25);
    }
}
