/***************************************************************************//**
 * Project: Colony
 *
 * \file    ChemicalSystemInitParam.h
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

#ifndef CHEMICALSYSTEMINITPARAM_H
#define CHEMICALSYSTEMINITPARAM_H

#include "../debug.h"
#include <blitz/array.h>
using blitz::Array;


/**
 * Initialize parameters structure for ChemicalSystem class.
 */
class ChemicalSystemInitParam
{
public:

  int nSpecies_; ///< Number of chemical species in the cells.
  int nChannels_; ///< Number of chemical reaction channels in the cells.
  Array<int,2> stoichMatrix_; ///< Stoichiometrix matrix defining the reactions in the cells.
  void (*computePropensitiesPtr_) (const Array<int,1>& x,
                                   Array<double,1>& propensities);
  /**< Pointer to the function that computes the propensities of the reaction channels in the cells.*/
  void (*computePropensitiesTimeDependentPtr_)
  (
    const Array<int,1>& x,
    const double volume,
    const double volume0,
    Array<double,1>& propensities
  );
  /**< Pointer to the function that computes the propensities of the reaction channels in the cells, time-dependent version.*/
};

#endif // CHEMICALSYSTEMINITPARAM_H
