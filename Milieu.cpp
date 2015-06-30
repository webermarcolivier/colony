/***************************************************************************//**
 * Project: Colony
 *
 * \file    Milieu.cpp
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

#include "Milieu.h"

//------------------------------------------------------------------------------

double Milieu::getTimeDependentVolume(const double time) const
{
  return milieuVolume_;
}

//------------------------------------------------------------------------------


