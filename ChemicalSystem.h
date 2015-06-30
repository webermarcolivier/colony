/***************************************************************************//**
 * Project: Colony
 *
 * \file    ChemicalSystem.h
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

#ifndef CHEMICALSYSTEM_H
#define CHEMICALSYSTEM_H

#include "debug.h"

// libraries header files
#include <blitz/array.h>
#include "Param/ChemicalSystemInitParam.h"

// namespaces
using blitz::Array;


/**
  * Chemical system of the cell (reactions).
  *
  * This class defines the chemical system of a cell. It defines all the reaction
  * channels through the stoichiometric matrix, and the function computing the
  * propensities of reaction channels.
  */
class ChemicalSystem
{
public:

//---- LIFECYCLE

  ChemicalSystem ();

  ~ChemicalSystem ();

  void initialize(const ChemicalSystemInitParam& p);


//---- OPERATORS

  /**
   * Assignment operator.
   */
  ChemicalSystem& operator = (const ChemicalSystem& c2);


public:

//---- DATA

// Rem: All data are PUBLIC.

  /**
   * Number of chemical species.
   */
  int nSpecies_;

  /**
   * Number of reaction channels.
   */
  int nChannels_;

  /**
   * Stoichiometric matrix of the chemical reactions.
   * Its dimension should be #nChannels x #nSpecies.
   */
  Array<int,2> stoichMatrix_;

  /**
   * Pointer to the function that computes the propensities of the different
   * reaction channels. The function is meant to be written by an external
   * program. Otherwise, it can be directly included in the ChemicalSystem
   * class.
   */
  void (*computePropensitiesPtr_)  (const Array<int,1>& x,
                                   Array<double,1>& propensities);

  /**
   * Pointer to the function that computes the propensities of the different
   * reaction channels. The function is meant to be written by an external
   * program. Otherwise, it can be directly included in the ChemicalSystem
   * class.
   * Time-dependent version. For example, for a variable volume with second-order
   * reactions.
   */
  void (*computePropensitiesTimeDependentPtr_)
  (
    const Array<int,1>& x,
    const double volume,
    const double volume0,
    Array<double,1>& propensities
  );


private:

  /**
   * Preventing copy constructor use.
   */
   ChemicalSystem(const ChemicalSystem& c2);

  /**
   * Copy method used in the assigment operator and the copy constructor.
   */
  void copy(const ChemicalSystem& c2);


};

#endif // CHEMICALSYSTEM_H
