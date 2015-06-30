/***************************************************************************//**
 * Project: Colony
 *
 * \file    IntegratorContext.cpp
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

#include "IntegratorContext.h"

//------------------------------------------------------------------------------

template<class IntegratorType>
IntegratorContext<IntegratorType>::IntegratorContext(Simulator* simulatorPtr)
  : integrator_(simulatorPtr)
{}

//------------------------------------------------------------------------------

template<class IntegratorType>
IntegratorContext<IntegratorType>::~IntegratorContext()
{}

//------------------------------------------------------------------------------

/**
 * To avoid linker errors, we have to define the template classes we will use.
 * For more information, see the C++ FAQ:\n
 * http://www.parashift.com/c++-faq-lite/templates.html#faq-35.15
 */
template class IntegratorContext<IntegratorGillespie>;
template class IntegratorContext<IntegratorGillespieModified>;
#ifdef USE_CHEMICAL_LANGEVIN
  template class IntegratorContext<IntegratorChemicalLangevin>;
#endif

//------------------------------------------------------------------------------
