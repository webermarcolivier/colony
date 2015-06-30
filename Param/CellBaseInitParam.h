/***************************************************************************//**
 * Project: Colony
 *
 * \file    CellBaseInitParam.h
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

#ifndef CELLBASEINITPARAM_H
#define CELLBASEINITPARAM_H

#include "../debug.h"
#include "StateInitParam.h"
#include "ChemicalSystemInitParam.h"



/**
 * Initialize parameters structure for CellBase class.
 */
class CellBaseInitParam
{
public:

  StateInitParam stateInitParam_;

  ChemicalSystemInitParam chemicalSystemInitParam_;
};


#endif // CELLBASEINITPARAM_H
