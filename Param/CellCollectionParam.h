/***************************************************************************//**
 * Project: Colony
 *
 * \file    CellCollectionParam.h
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

#ifndef CELLCOLLECTIONPARAM_H
#define CELLCOLLECTIONPARAM_H

#include "../debug.h"

// standard C++ header files

// libraries header files
#include <blitz/array.h>

// user header files
#include "CellBaseInitParam.h"
#include "CellInitParam.h"

// namespaces



class CellCollectionParam
{
public:

  int nCells_; ///< Number of cells in the cell collection.

  double totalVolume_; ///< Total volume of the system (milieu + cells).

  CellInitParam cellInitParam_; ///< Parameters for initializing the cells.

  blitz::Array<int,2> x0Cells; ///< Initial conditions for all cells (number of molecules of each species in every cell).

  CellBaseInitParam milieuInitParam_; ///< Parameters for initializing the milieu.

  bool constantCellDensity_; ///< Keep the cell density constant by automatically killing cells when duplicating.

  bool x0CellHeterogenous_;

  Array<int,2> x0CellHeterogenousTable_;

  Array<double,1> initialPhaseTable_;

};

#endif // CELLCOLLECTIONPARAM_H
