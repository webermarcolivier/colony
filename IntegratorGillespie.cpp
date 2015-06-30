/***************************************************************************//**
 * Project: Colony
 *
 * \file    IntegratorGillespie.cpp
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

#include "IntegratorGillespie.h"

//------------------------------------------------------------------------------

IntegratorGillespie::IntegratorGillespie(Simulator* simulatorPtr)
  : simulatorPtr_(simulatorPtr)
{}

//------------------------------------------------------------------------------

IntegratorGillespie::~IntegratorGillespie()
{}

//------------------------------------------------------------------------------

void IntegratorGillespie::integrateOneStep()
{
  ///- Instantiate a local reference to the array of propensities.
  Array<double,1> propensities
  (
    simulatorPtr_->cellCollection_.computeAndGetPropensities()
  );

  ///- Compute the total sum of all channels propensities, 'atot'.
  double atot = simulatorPtr_->cellCollection_.getSumPropensities();

  #ifdef BUILD_DEBUG
  // We should check that the optimized computation of the sum of propensities
  // does not diverge from the normal sum.
  // Accumulation of addition and substraction -> accumulation of numerical
  // imprecision.
  double atot2 = sum( propensities );
  if (abs(atot - atot2) > 1e-09)
  {
    cout << "WARNING: class IntegratorGillespie, function integrateOneStep(),"
        " optimized computation of the sum of propensities diverges from the normal sum."
        << endl;
  }
  #endif // BUILD_DEBUG

  ///- Pick up 2 uniform random numbers in [0,1), 'r1' and 'r2'.
  double r1 = RandomNumberGenerator::getUniform();
  double r2 = RandomNumberGenerator::getUniform();

  ///- Compute the time till next reaction, 'tau'.
  double tau;
  if (atot != 0.0) {
    tau = log(1.0/r1)/atot;
  } else {
    tau = numeric_limits<double>::infinity();
  }

  ///- Find the channel of the next reaction, 'mu'.
  ///  Choose channel between all the reactions available according to the weights
  ///  given by the propensities.
  int channel = 0;
  double asum = propensities(0);
  int nChannels = propensities.extent(firstDim);

  while (asum <= atot*r2 && channel < nChannels-1)
    /* The condition channel < nChannels-1 is not necessary but is there to
       prevent numerical problems with the sum of propensities. For example,
       if the sum atot2 is not slightly inferior to the real sum, then we could
       have a channel number greater than the real number of channels. */
  {
    channel++;
    asum += propensities(channel);
  }

  ///- Get the time of the following division event in the whole cell collection.
  ///  The next division event will happen at time 't1'.
  double timeNextDivision =
    simulatorPtr_->cellCollection_.getTimeNextDivision();

  ///- If the time of next reaction is smaller than the time of next division
  ///  event 't1', then
  ///  - Apply reaction 'mu' in cell 'k' and update simulation time.
  ///- If not, update simulation time to the time of next division
  ///  event 't1' and divide the cell(s).
  if (simulatorPtr_->getTime() + tau <= timeNextDivision)
  {
    simulatorPtr_->increaseTime(tau);

    // Update state of the system by the stoichiometric vector of the
    // corresponding channel.
    simulatorPtr_->cellCollection_.applyReaction(channel);

  } else {
    // Jump to the next cell cycle.
    simulatorPtr_->setTime(timeNextDivision);

    // Duplicate cell(s).
    simulatorPtr_->cellCollection_.applyNextDivisionEvent();
  }

}



