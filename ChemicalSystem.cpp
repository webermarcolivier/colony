/***************************************************************************//**
 * Project: Colony
 *
 * \file    ChemicalSystem.cpp
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

#include "ChemicalSystem.h"

//------------------------------------------------------------------------------

ChemicalSystem::ChemicalSystem()
  : nSpecies_(0),
    nChannels_(0),
    stoichMatrix_(),  //Rem: the Array is 0-initialized, its data cannot be accessed
    computePropensitiesPtr_(0),
    computePropensitiesTimeDependentPtr_(0)
{}

//------------------------------------------------------------------------------

ChemicalSystem::~ChemicalSystem()
{
  stoichMatrix_.free();
}

//------------------------------------------------------------------------------

void ChemicalSystem::initialize(const ChemicalSystemInitParam& p)
{
  nSpecies_ = p.nSpecies_;
  nChannels_ = p.nChannels_;

  stoichMatrix_.resize( p.stoichMatrix_.rows(), p.stoichMatrix_.cols() );
  stoichMatrix_ = p.stoichMatrix_;

  #ifndef TIME_DEPENDENT_PROPENSITIES
    computePropensitiesPtr_ = p.computePropensitiesPtr_;
  #else
    computePropensitiesTimeDependentPtr_ = p.computePropensitiesTimeDependentPtr_;
  #endif // TIME_DEPENDENT_PROPENSITIES
}

//------------------------------------------------------------------------------

void ChemicalSystem::copy(const ChemicalSystem& c2)
{
  nSpecies_ = c2.nSpecies_;
  nChannels_ = c2.nChannels_;

  stoichMatrix_.resize( c2.stoichMatrix_.rows(), c2.stoichMatrix_.cols() );
  stoichMatrix_ = c2.stoichMatrix_;

  computePropensitiesPtr_ = c2.computePropensitiesPtr_;
  computePropensitiesTimeDependentPtr_ = c2.computePropensitiesTimeDependentPtr_;
}

//------------------------------------------------------------------------------

ChemicalSystem& ChemicalSystem::operator = (const ChemicalSystem& c2)
{
  // Test for self-assignment
  if (&c2 == this)
  {
    // Do nothing and just return myself.
    return *this;
  }

  copy(c2);

  return *this;
}

//------------------------------------------------------------------------------



