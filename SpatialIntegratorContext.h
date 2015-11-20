/***************************************************************************//**
 * Project: Colony
 *
 * \file    SpatialIntegratorContext.h
 * \author  Marc Weber\n
 *          The SiMBioSys group (CosmoLab)\n
 *          Parc Científic de Barcelona\n
 *          Barcelona, Spain.\n
 *          http://www.thesimbiosys.com
 * \version 1.0
 * \date    05/2011
 *
 *          Copyright 2009 by Marc Weber
 ******************************************************************************/

#ifndef SPATIALINTEGRATORCONTEXT_H
#define SPATIALINTEGRATORCONTEXT_H

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
  * SpatialIntegratorContext class is a template class with template type defining the
  * type of spatial integrator.
  * It contains an SpatialIntegratorType object and uses its integrateOneStep() function.
  */
template<class SpatialIntegratorType>
class SpatialIntegratorContext
{
public:

  /**
   * Initiate the SpatialIntegratorContext with the pointer to the Simulator object
   * as argument.
   */
  SpatialIntegratorContext(Simulator* simulatorPtr);

  ~SpatialIntegratorContext();

  SpatialIntegratorType integrator_;

  void integrate(double timeStep)
  {
    integrator_.integrate(timeStep);
  }


  void setIsEquilibrationSteps(bool isEquilibrationSteps)
  {
    integrator_.setIsEquilibrationSteps(isEquilibrationSteps);
  }
};

#include "Simulator.h"

#endif // SPATIALINTEGRATORCONTEXT_H
