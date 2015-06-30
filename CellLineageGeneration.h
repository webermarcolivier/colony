/***************************************************************************//**
 * Project: Colony
 *
 * \file    CellLineageGeneration.h
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

#ifndef CELLLINEAGEGENERATION_H
#define CELLLINEAGEGENERATION_H

#include "debug.h"

// standard C++ header files
#include <iostream>

// libraries header files
#include <blitz/array.h>

// user header files

// namespaces
using blitz::Array;
using std::cout;
using std::ostream;
using std::endl;



/**
 * %Cell lineage event defining the genealogy of cells after a change in cell number.
 */
class CellLineageGeneration
{
public:

  CellLineageGeneration();
  ~CellLineageGeneration();
  CellLineageGeneration(const CellLineageGeneration& t);

  void initialize
  (
    double time,
    int nCells,
    int nCellsDelta,
    Array<int,1> motherCellsIndices
  );

  /**
   * Time of the cell lineage event (division(s), cell deletion).
   */
  double time_;

  /**
   * Change in the number of cells.
   */
  int nCellsDelta_;

  /**
   * Number of cells *after* the lineage event.
   */
  int nCells_;

  /**
   * Vector containing a list of the cells indices transformation.
   * Element i of the vector contains the index of the parent cell for cell i.
   * For a division, the vector will contains two elements with index i0,
   * where i0 is the index (in the array before the division) of the mother cell.
   *
   * Example:
   * (0, 1, 2), division of cell 1 => (0, 1, 2, 1).
   *
   * For a deletion, the vector contains the indices of the rearranged cells
   * in the new array with one element less.
   *
   * Example:
   * (0, 1, 2, 3), deletion of cell 1 => (0, 2, 3).
   */
  Array<int,1> motherCellsIndices_;
};




#endif // CELLLINEAGEGENERATION_H
