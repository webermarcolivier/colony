/***************************************************************************//**
 * Project: Colony
 *
 * \file    IntegratorChemicalLangevin.cpp
 * \author  Marc Weber\n
 *          The SiMBioSys group (CosmoLab)\n
 *          Parc Científic de Barcelona\n
 *          Barcelona, Spain.\n
 *          http://www.thesimbiosys.com
 * \version 1.0
 * \date    03/2012
 *
 *          Copyright 2012 by Marc Weber
 ******************************************************************************/

#include "debug.h"

#ifdef USE_CHEMICAL_LANGEVIN

#include "IntegratorChemicalLangevin.h"

//------------------------------------------------------------------------------

IntegratorChemicalLangevin::IntegratorChemicalLangevin(Simulator* simulatorPtr)
  : simulatorPtr_(simulatorPtr)
{}

//------------------------------------------------------------------------------

IntegratorChemicalLangevin::~IntegratorChemicalLangevin()
{}

//------------------------------------------------------------------------------

void IntegratorChemicalLangevin::integrateOneStep(double timeStep)
{
  chemicalLangevinComputeIncrement_.chemicalLangevinComputeIncrement(
    simulatorPtr_->cellCollection_.getGlobalXConcCells(),
    simulatorPtr_->cellCollection_.getMilieu().getXConc(),
    simulatorPtr_->cellCollection_.getNCells(),
    simulatorPtr_->cellCollection_[0].getVolume0(),
    simulatorPtr_->cellCollection_.getMilieu().getVolume0(),
    timeStep
  );

  simulatorPtr_->increaseTime(timeStep);


  simulatorPtr_->proportionNegativeSteps_ = chemicalLangevinComputeIncrement_.proportionNegativeSteps;

  // In the case of growing and dividing cells:
  // Check the time of the following division event and adapt time step,
  // divide cells, etc.

}

#endif //USE_CHEMICAL_LANGEVIN
