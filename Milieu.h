/***************************************************************************//**
 * Project: Colony
 *
 * \file    Milieu.h
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

#ifndef MILIEU_H
#define MILIEU_H

#include "compilation_options.h"

// libraries header files
#include <blitz/array.h>

// user header files
#include "CellBase.h"

// namespaces
using blitz::Array;
using blitz::Range;
using std::cout;
using std::ostream;
using std::endl;



/**
  * A special cell that represents the environment in which the
  * colony of cells live. It can contain chemical species and reactions, like
  * a normal CellBase object and cannot divide.
  */

class Milieu : public CellBase
{
public:

  /**
   * We reimplement the method in order to return a value of the milieu volume
   * independent of time and dependent only on the number of cells (volume0).
   * This is an approximation since we should take into account the volume of
   * all the other cells in order to calculate the volume of the milieu. However,
   * as the volume of the milieu is much larger than the volume of one cell,
   * there is not much difference.
   * UPDATE 2013.10.17: change this method in order to take into account the volume of
   * ALL the cells and to compute exactly the volume of the milieu. Managed by CellCollection
   * (this is easier since we have to know the volume of ALL the cells).
   * We use a variable milieuVolume_ that we will only update at each time slice,
   * trade-off between accuracy and computation time. Also,
   * volume0 is now constant for the milieu and refers to the inital volume.
   */
  double getTimeDependentVolume(const double time) const;

  double milieuVolume_;

};

#endif // MILIEU_H
