/***************************************************************************//**
 * Project: Colony
 *
 * \file    IntegratorChemicalLangevin.h
 * \author  Marc Weber\n
 *          The SiMBioSys group (CosmoLab)\n
 *          Parc Científic de Barcelona\n
 *          Barcelona, Spain.\n
 *          http://www.thesimbiosys.com
 * \version 1.0
 * \date    03/2012
 *
 *          Copyright 2012 by Marc Weber
 ******************************************************************************/

#ifndef INTEGRATORCHEMICALLANGEVIN_H
#define INTEGRATORCHEMICALLANGEVIN_H

#include "debug.h"

// standard C++ header files
#include <iostream>
#include <limits> ///< numerical limits for each of the fundamental types

// libraries header files
#include <blitz/array.h>

// user header files
#include "RandomNumberGenerator.h"
#include "chemicalLangevinComputeIncrement.h"

/**
 * Use forward declaration of class Simulator. The header file "Simulator.h"
 * is included after the declaration of the IntegratorGillespie class.
 * This is necessary since the IntegratorGillespie class contains a pointer to
 * the Simulator class.
 */
class Simulator;


// namespaces
using blitz::Array;
using std::cout;
using std::ostream;
using std::endl;
using std::numeric_limits;



/**
  * Algorithm of Euler for the chemical Langevin equation.
  */
class IntegratorChemicalLangevin
{
public:

  /**
   * Initiate the IntegratorContext with the pointer to the Simulator object
   * as argument.
   */
  IntegratorChemicalLangevin(Simulator* simulatorPtr);

  ~IntegratorChemicalLangevin();

  /**
   * This method is not used in this integrator class.
   */
  void integrateOneStep()
  {}

  /**
   * Integrate one step in the chemical Langevin algorithm.
   */
  void integrateOneStep(double timeStep);

  void setFixedNCells(int nCells)
  {
    chemicalLangevinComputeIncrement_.setFixedNCells(nCells);
  }


private:

  Simulator* simulatorPtr_;

  ChemicalLangevinComputeIncrement chemicalLangevinComputeIncrement_;

  /**
   * Pointer to the function that computes one step of the Langevin
   * algorithm. The function is meant to be written by an external
   * program. Otherwise, it can be directly included in this
   * class.
   */
//  void (*chemicalLangevinComputeIncrementPtr_)
//  (
//    Array<double,1>& xConcCells,
//    Array<double,1>& xConcMilieu,
//    const int nCells,
//    const double volumeCell,
//    const double volumeExt,
//    const double chemicalLangevinTimeStep
//  );

};

//------------------------------------------------------------------------------

#include "Simulator.h"

#endif // INTEGRATORCHEMICALLANGEVIN_H
