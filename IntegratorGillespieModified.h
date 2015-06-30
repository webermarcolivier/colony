/***************************************************************************//**
 * Project: Colony
 *
 * \file    IntegratorGillespieModified.h
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

#ifndef INTEGRATORGILLESPIEMODIFIED_H
#define INTEGRATORGILLESPIEMODIFIED_H

#include "debug.h"

// standard C++ header files
#include <iostream>
#include <limits> ///< numerical limits for each of the fundamental types

// libraries header files

// user header files
#include "RandomNumberGenerator.h"

/**
 * Use forward declaration of class Simulator. The header file "Simulator.h"
 * is included after the declaration of the class.
 * This is necessary since the class contains a pointer to
 * the Simulator class.
 */
class Simulator;


// namespaces
using std::cout;
using std::ostream;
using std::endl;
using std::numeric_limits;



/**
  * Algorithm of Gillespie modified for time-dependent diffusion process.
  *
  * The algorithm is using the approximate method that consists of taking the
  * time till next reaction (\f$ \tau \f$) distribution of the normal Gillespie
  * algorithm to compute \f$ \tau \f$ and picking the reaction channel based on
  * the recalculated propensities at time \f$ t+ \tau \f$.
  */
class IntegratorGillespieModified
{
public:

  /**
   * Initiate the IntegratorContext with the pointer to the Simulator object
   * as argument.
   */
  IntegratorGillespieModified(Simulator* simulatorPtr);

  ~IntegratorGillespieModified();

  /**
   * Integrate one step in the Gillespie algorithm.
   */
  void integrateOneStep();

  /**
   * This method is not used in this integrator class.
   */
  void integrateOneStep(double timeStep)
  {}

private:

  Simulator* simulatorPtr_;

  double nStep_;

};

//------------------------------------------------------------------------------

#include "Simulator.h"



#endif // INTEGRATORGILLESPIEMODIFIED_H
