/***************************************************************************//**
 * Project: Colony
 *
 * \file    IntegratorGillespie.h
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

#ifndef INTEGRATORGILLESPIE_H
#define INTEGRATORGILLESPIE_H

#include "compilation_options.h"

// standard C++ header files
#include <iostream>
#include <limits> ///< numerical limits for each of the fundamental types

// libraries header files

// user header files
#include "RandomNumberGenerator.h"

/**
 * Use forward declaration of class Simulator. The header file "Simulator.h"
 * is included after the declaration of the IntegratorGillespie class.
 * This is necessary since the IntegratorGillespie class contains a pointer to
 * the Simulator class.
 */
class Simulator;


// namespaces
using std::cout;
using std::ostream;
using std::endl;
using std::numeric_limits;



/**
  * Algorithm of Gillespie.
  */
class IntegratorGillespie
{
public:

  /**
   * Initiate the IntegratorContext with the pointer to the Simulator object
   * as argument.
   */
  IntegratorGillespie(Simulator* simulatorPtr);

  ~IntegratorGillespie();

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

};

//------------------------------------------------------------------------------

#include "Simulator.h"

#endif // INTEGRATORGILLESPIE_H
