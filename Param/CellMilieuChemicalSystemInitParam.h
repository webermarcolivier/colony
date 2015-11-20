/***************************************************************************//**
 * Project: Colony
 *
 * \file    CellMilieuChemicalSystemInitParam.h
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

#ifndef CELLMILIEUCHEMICALSYSTEMINITPARAM_H
#define CELLMILIEUCHEMICALSYSTEMINITPARAM_H

#include "../compilation_options.h"
#include <blitz/array.h>
using blitz::Array;
class Milieu;
class Cell;



/**
 * Initialize parameters structure for CellMilieuChemicalSystem class.
 */
class CellMilieuChemicalSystemInitParam
{
public:

  Milieu* milieuPtr_;
  Cell* cellPtr_;
  int nChannelsCellMilieu_; ///< Number of cell-milieu reaction channels.
  Array<int,2> stoichMatrixCellMilieu_; ///< Stoichiometrix matrix defining the cell-milieu reactions.
  /**
   * Pointer to the function that computes the propensities of the cell-milieu
   * reaction channels.
   * @see CellMilieuChemicalSystem::computePropensitiesInterCellularPtr_
   */
  void (*computePropensitiesInterCellularPtr_)
    (const Array<int,1>&,
     const Array<int,1>&,
     const double,
     const double,
     Array<double,1>&);
  /**
   * Pointer to the function that computes the time-dependent propensities of
   * the cell-milieu reaction channels.
   * @see CellMilieuChemicalSystem::computeTimeDependentPropensitiesInterCellularPtr_
   */
  void (*computeTimeDependentPropensitiesInterCellularPtr_)
  (
    const Array<int,1>& x1,
    const Array<int,1>& x2,
    const double volume,
    const double volume0,
    const double volumeExt,
    Array<double,1>& a
  );

};

#endif // CELLMILIEUCHEMICALSYSTEMINITPARAM_H
