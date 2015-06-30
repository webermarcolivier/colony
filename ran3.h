#pragma once

#define MBIG 1000000000 /* According to Knuth, any large MBIG, and */
#define MSEED 161803398 /* any smaller (but still large MSEED can  */
                        /* be substituted for the above values.    */       
#define MZ 0           
#define FAC (1.0/MBIG)

float ran3(long *idnum);
double unif_rand_dbl(long *idum);
double box_muller(double m, double s, long *idnum);
double expdev(double tau, long *idnum);

//long llavor;
