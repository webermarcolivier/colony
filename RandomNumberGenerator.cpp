/***************************************************************************//**
 * Project: Colony
 *
 * \file    RandomNumberGenerator.cpp
 * \author  Marc Weber\n
 *          The SiMBioSys group (CosmoLab)\n
 *          Parc Científic de Barcelona\n
 *          Barcelona, Spain.\n
 *          http://www.thesimbiosys.com
 * \version 1.0
 * \date    11/2009
 *
 *          Copyright 2009 by Marc Weber
 ******************************************************************************/

#include "RandomNumberGenerator.h"

//------------------------------------------------------------------------------

gsl_rng* RandomNumberGenerator::rngPtr_ = 0;
long RandomNumberGenerator::idnum_ = -1;
int RandomNumberGenerator::nLuxIExp_ = 0;

//------------------------------------------------------------------------------

RandomNumberGenerator::RandomNumberGenerator()
{}

//------------------------------------------------------------------------------

RandomNumberGenerator::~RandomNumberGenerator()
{}

//------------------------------------------------------------------------------

void RandomNumberGenerator::initialize()
{
  if (rngPtr_ == 0)
  {
    // Initialize random number generator of type ranlux.
    rngPtr_ = gsl_rng_alloc (gsl_rng_ranlxd1);
  }

  // Initialize r.n.g. seed. with clock generated seed
  unsigned long int rng_seed = (unsigned long int)time(0);
  gsl_rng_set(rngPtr_, rng_seed);
}

//------------------------------------------------------------------------------

void RandomNumberGenerator::initialize(int seed)
{
  if (rngPtr_ == 0)
  {
    // Initialize random number generator of type ranlux.
    rngPtr_ = gsl_rng_alloc (gsl_rng_ranlxd1);
    //rngPtr_ = gsl_rng_alloc (gsl_rng_mt19937);
  }

  // Initialize r.n.g. seed. with user defined seed
  gsl_rng_set(rngPtr_, seed);
  //gsl_rng_set(rngPtr_, 0);

  idnum_ = seed;
  unif_rand_dbl(&idnum_);
}

//------------------------------------------------------------------------------

double RandomNumberGenerator::getCellCycle(double tdet, double gamma)
{
  double tc;
  double tsto;

  tsto = gsl_ran_exponential(rngPtr_, tdet);
  tc   = gamma*tdet + (1.0 - gamma)*tsto;

  return tc;
}

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// RAN3.C implementation--------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

 /*****************************************************************************
* Generador de numeros aleatorios ran3                                       *
* Tomado de Numerical recipes in C de Press et al.                           *
*                                                                            *
* Esta funcion genera numeros aleatorios con distribucion uniforme en el     *
* rango [0.0, 1.0].  La rutina puede convertirse completamente a artimetica  *
* flotante declarando mj, mk, ma[] como float, y definiendo mbig y mseed     *
* como 4000000 y 1618033, respectivamente.                                   *
*                                                                            *
* Manuel Valenzuela (23 enero 1996)                                          *
 *****************************************************************************/




float ran3(long *idnum)
/********************************************************************
* Returns a uniform random deviate between 0.0 and 1.0. Set idnum   *
* to any negative value to initialize or reinitialize the sequence. *
********************************************************************/
{
	static int inext, inextp;
	static long ma[56];   /* The value 56 (range ma[1..55]) is special */
	/* and should not be modified; see Knuth.    */
	static int iff=0;
	long mj, mk;
	int i, ii, k;

	if (*idnum<0 || iff==0) {    /* Initialization */
		iff = 1;
		/* Initialize ma[55] using the seed idnum and the large number MSEED */
		mj = MSEED - (*idnum<0 ? -*idnum : *idnum);
		mj %= MBIG;
		ma[55] = mj;
		mk = 1;
		/* Now initizalize the rest of the table, in a slightly       */
		/* random order, with numbers that are not especially random. */
		for (i=1; i<=54; i++) {
			ii = (21*i) % 55;
			ma[ii] = mk;
			mk = mj - mk;
			if (mk < MZ) mk += MBIG;
			mj = ma[ii];
		}
		/* We randomize them by "warming up the generator." */
		for (k=1; k<=4; k++)
			for (i=1; i<=55; i++) {
			ma[i] -= ma[1+(i+30) % 55];
			if (ma[i]<MZ) ma[i] += MBIG;
			}
			inext = 0;     /* Prepare indices for our first generated number. */
			inextp = 31;   /* The constant 31 is special; see Knuth */
			*idnum = 1;
	}
	/* Here is where we start, except on initialization */
	if (++inext==56) inext = 1;   /* Initizalize inext and inextp, wrapping    */
	if (++inextp==56) inextp = 1; /* arround 56 to 1.                          */
	mj = ma[inext] - ma[inextp];  /* Generate new random number substractively */
	if (mj<MZ) mj +=MBIG;         /* Make sure that it is in range.            */
	ma[inext] = mj;               /* Store it                                  */


	return mj*FAC;                /* and output the derived uniform deviate.   */
}

/************************************************************************/
/* unif_rand_dbl() returns a uniform random variate between [0.0, 1.0]. */
/*      (based on NumRec routine rand3(), based itself on this-or-that  */
/*      from Knuth.  Not a linear congruence generator.                 */
/* To initialize/reinitialize, pass it a negative long int; it has a    */
/*      memory, so passing it the same initializer multiple times       */
/*      during a run of the program will produce different values.      */
/************************************************************************/
double unif_rand_dbl(long *idnum)
{
	double highorder = (double) ran3(idnum);
	double loworder = (double) ran3(idnum);

	return highorder + loworder*FAC;
}

double box_muller(double m, double s, long *idnum)	/* normal random variate generator */
{				        /* mean m, standard deviation s */
	double x1, x2, w, y1;
	static double y2;
	static int use_last = 0;

	if (use_last)		        /* use value from previous call */
	{
		y1 = y2;
		use_last = 0;
	}
	else
	{
		do {
			x1 = 2.0 * unif_rand_dbl(idnum) - 1.0;
			x2 = 2.0 * unif_rand_dbl(idnum) - 1.0;
			w = x1 * x1 + x2 * x2;
		} while ( w >= 1.0 );

		w = sqrt( (-2.0 * log( w ) ) / w );
		y1 = x1 * w;
		y2 = x2 * w;
		use_last = 1;
	}

	//printf("%f",s);
	return( m + y1 * s );

}

double expdev(double tau, long *idnum){

	double dum;

	do

	    dum=unif_rand_dbl(idnum);

	while (dum == 0.0);

	//fprintf(stderr,"\ndum=%lf joder=%lf tau=%lf\n",dum,joder,tau);

	return -tau*log(dum);

}

//------------------------------------------------------------------------------




