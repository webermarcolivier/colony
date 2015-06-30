/***************************************************************************//**
 * Project: Colony
 *
 * \file    SpatialIntegratorContext.cpp
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

#include "SpatialIntegratorContext.h"

//------------------------------------------------------------------------------

template<class SpatialIntegratorType>
SpatialIntegratorContext<SpatialIntegratorType>::SpatialIntegratorContext(Simulator* simulatorPtr)
  : integrator_(simulatorPtr)
{}

//------------------------------------------------------------------------------

template<class SpatialIntegratorType>
SpatialIntegratorContext<SpatialIntegratorType>::~SpatialIntegratorContext()
{}

//------------------------------------------------------------------------------

/**
 * To avoid linker errors, we have to define the template classes we will use.
 * For more information, see the C++ FAQ:\n
 * http://www.parashift.com/c++-faq-lite/templates.html#faq-35.15
 */
template class SpatialIntegratorContext<SpatialIntegratorODE>;

//------------------------------------------------------------------------------

