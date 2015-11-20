/***************************************************************************//**
 * Project: Colony
 *
 * \file    CellInitParam.h
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

#ifndef CELLINITPARAM_H
#define CELLINITPARAM_H

#include "../compilation_options.h"
#include "GraphicsCellCompositeParam.h"
#include "CellBaseInitParam.h"
#include "CellMilieuChemicalSystemInitParam.h"

/**
 * Initialize parameters structure for Cell class.
 */
class CellInitParam
{
public:

  CellBaseInitParam cellBaseInitParam_;

  double cellCycleDeterministic_; ///< Deterministic value of the duration of the cell cycle.
  double cellCycleGamma_; ///< Stochastic to deterministic weight parameter for the cell cycle distribution.
  Array<int,2> dnaSpecies_; ///< Matrix defining the genetic composition of the species. @see Cell#dnaSpecies_

  double cellCyclephase_; ///< Initial phase of the cell cycle, values in [0,1].

  CellMilieuChemicalSystemInitParam cellMilieuChemicalSystemInitParam_;
  // Do not forget to initiate the cellPtr_ member!

  GraphicsCellCompositeParam graphicsCellParam_;
};


#endif // CELLINITPARAM_H
