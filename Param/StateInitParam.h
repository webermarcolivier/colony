/***************************************************************************//**
 * Project: Colony
 *
 * \file    StateInitParam.h
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

#ifndef STATEINITPARAM_H
#define STATEINITPARAM_H

#include "../compilation_options.h"
#include <vector>
#include <string>
#include "blitz/array.h"
using blitz::Array;
using std::vector;
using std::string;


/**
 * Initialize parameters structure for State class.
 */
class StateInitParam
{
public:

  Array<int,1> x_; ///< Initial state of the cells (number of molecules of chemical species).
  vector<string> listSpeciesName_; ///< List of the species names in the cells.
  int propensitiesSize_;
  double cellVolume0_; ///< Initial volume of the cells.
  double cellLength0_;
  double cellHeight0_;
};

#endif // STATEINITPARAM_H
