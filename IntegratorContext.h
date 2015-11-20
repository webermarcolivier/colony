/***************************************************************************//**
 * Project: Colony
 *
 * \file    IntegratorContext.h
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

#ifndef INTEGRATORCONTEXT_H
#define INTEGRATORCONTEXT_H

#include "compilation_options.h"

// standard C++ header files
#include <iostream>

// libraries header files

// user header files

// namespaces
using std::cout;
using std::ostream;
using std::endl;

/**
 * Use forward declaration of class Simulator. The header file "Simulator.h"
 * is included after the declaration of the IntegratorContext class.
 * This is necessary since the IntegratorContext class contains a pointer to
 * the Simulator class.
 */
class Simulator;



/**
  * Template class that defines the context for the integrator. The
  * IntegratorContext class is a template class with template type defining the
  * type of integrator (Gillespie, modified Gillespie, etc).
  * It contains an IntegratorType object and uses its integrateOneStep() function.
  */
template<class IntegratorType>
class IntegratorContext
{
public:

  /**
   * Initiate the IntegratorContext with the pointer to the Simulator object
   * as argument.
   */
  IntegratorContext(Simulator* simulatorPtr);

  ~IntegratorContext();


  void integrateOneStep()
  {
    integrator_.integrateOneStep();
  }

  void integrateOneStep(double timeStep)
  {
    integrator_.integrateOneStep(timeStep);
  }

  IntegratorType integrator_;


};

#include "Simulator.h"

#endif // INTEGRATORCONTEXT_H
