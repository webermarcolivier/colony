/***************************************************************************//**
 * Project: Colony
 *
 * \file    RandomNumberGenerator.h
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

#ifndef RANDOMNUMBERGENERATOR_H
#define RANDOMNUMBERGENERATOR_H

#include "debug.h"

// standard C++ header files
#include <time.h>

// libraries header files
#include <gsl/gsl_rng.h> ///< random number generator from the GSL library.
#include <gsl/gsl_randist.h> ///< random number distribution from the GSL library.

// user header files

// namespaces

//------------------------------------------------------------------------------
// RAN3.H definitions-----------------------------------------------------------
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
//------------------------------------------------------------------------------


/**
 * Random number generator from the GSL library.
 */
class RandomNumberGenerator
{
public:

//---- LIFECYCLE

  RandomNumberGenerator();

  ~RandomNumberGenerator();

  /**
   * Initiate the random number generator with a clock-cycle based seed.
   */
  static void initialize();

  /**
   * Initiate the random number generator with user defined seed.
   */
  static void initialize(int seed);


//---- OPERATORS


//---- ACCESS

  /**
   * Returns a double in [0,1) with uniform probability.
   */
  static double getUniform();

  /**
   * Returns a random number with normal distribution N(0,1).
   */
  static double getNormal();

  /**
   * Returns a random integer from the binomial distribution,
   * the number of successes in \c n independent trials with probability \c p.
   */
  static int getBinomial(double p, unsigned int n);

  /**
   * Returns a random value for the cell cycle \f$t_c\f$.
   * It is given by:
   * \f[ t_c = \gamma t_{det} + (1-\gamma)t_{sto} \f]
   * where \f$t_{det}\f$ and \f$t_{sto}\f$ are the deterministic and stochastic
   * components of the cell cycle respectively, and \f$\gamma\in [0,1]\f$ is a
   * parameter that weight their relative importance. Note that we acknowledge
   * the fact that, prior to division, cells must mature and grow and consequently
   * the cell cycle has a minimum duration \f$\gamma t_{det}\f$. The stochastic
   * part, \f$t_{sto}\f$, is supposed to be exponentially distributed and its
   * probability density reads,
   * \f[ \rho(t_{sto}) = \frac{1}{t_{det}} e^{-\frac{t_{sto}}{t_{det}}} \f]
   * Note that:
   * \f[ <\tau> = t_{det} \f]
   * \f[ \sigma_{\tau} = (1-\gamma) t_{det} \f]
   */
  static double getCellCycle(double tdet, double gamma);


//---- INQUIRY


//---- OPERATIONS


private:

//---- DATA

  /**
   * Pointer to the random number generator
   */
  static gsl_rng* rngPtr_;

  static long idnum_;

public:
  static int nLuxIExp_;
};

//------------------------------------------------------------------------------

inline double RandomNumberGenerator::getUniform()
  //{return unif_rand_dbl(&idnum_);}
  {return gsl_rng_uniform_pos(rngPtr_);}

inline double RandomNumberGenerator::getNormal()
  {return gsl_ran_ugaussian(rngPtr_);}

inline int RandomNumberGenerator::getBinomial(double p, unsigned int n)
  {return (int) gsl_ran_binomial (rngPtr_, p, n);}

//------------------------------------------------------------------------------

#endif // RANDOMNUMBERGENERATOR_H
