/***************************************************************************//**
 * Project: Colony
 *
 * \file    IntegratorGillespieModified.cpp
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

#include "IntegratorGillespieModified.h"

//------------------------------------------------------------------------------

IntegratorGillespieModified::IntegratorGillespieModified(Simulator* simulatorPtr)
  : simulatorPtr_(simulatorPtr),
    nStep_(0.0)
{}

//------------------------------------------------------------------------------

IntegratorGillespieModified::~IntegratorGillespieModified()
{}

//------------------------------------------------------------------------------

void IntegratorGillespieModified::integrateOneStep()
{
  nStep_++;

  ///- Instantiate a local reference to the array of propensities.
  ///- Compute all propensities.
  Array<double,1> propensities
  (
    simulatorPtr_->
      cellCollection_.computeAndGetPropensitiesModified(simulatorPtr_->getTime())
  );

  ///- Compute the total sum of all channels propensities at time t.
  double atot = simulatorPtr_->cellCollection_.getSumPropensities();

  #ifdef SUM_OF_PROPENSITIES_CHECK
  double errorTolerance = 1e-13;
  double atottrue = sum ( propensities );
  if (fabs(atot/atottrue - 1.0) > errorTolerance)
  {
    cout << "WARNING: class IntegratorGillespieModified, function integrateOneStep(),"
            " the optimized sum of propensities differs from the"
            " true calculated sum of propensities." << endl;
    cout << "Optimized sum of propensities = " << scientific << setprecision(15) << atot
         << ", Correct sum of propensities = " << scientific << setprecision(15) << atottrue
         << ", nStep = " << nStep_ << " time = " << simulatorPtr_->getTime() << endl;
  }
  #endif //SUM_OF_PROPENSITIES_CHECK


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

  ///- Recompute the propensities of time-dependent channels at time t + tau.
// TODO (mweber#1#): To implement: This function just re-calculate all the propensities, because there are time-dependent channels potentially everywhere. Should implement a function in Mathematica that calculate only the time-dep. channels inside a cell.
//  simulatorPtr_->cellCollection_.computeTimeDependentPropensities
//                                               (simulatorPtr_->getTime() + tau);
  ///- This is just to costly to do. Let the propensities like they are.


  ///- Compute the total sum of all channels propensities at time t + tau.
  double atot2 = simulatorPtr_->cellCollection_.getSumPropensities();

  #ifdef SUM_OF_PROPENSITIES_CHECK
  double atot2true = sum ( propensities );
  if (fabs(atot2/atot2true - 1.0) > errorTolerance)
  {
    cout << "WARNING: class IntegratorGillespieModified, function integrateOneStep(),"
            " the optimized sum of propensities differs from the"
            " true calculated sum of propensities." << endl;
    cout << "Optimized sum of propensities = " << scientific << setprecision(15) << atot
         << ", Correct sum of propensities = " << scientific << setprecision(15) << atottrue
         << ", nStep = " << nStep_ << " time = " << simulatorPtr_->getTime() << endl;
  }
  #endif //SUM_OF_PROPENSITIES_CHECK

  ///- Find the channel of the next reaction, 'mu'.
  ///  Choose channel between all the reactions available according to the weights
  ///  given by the propensities.
  int channel = 0;
  int nChannels = propensities.extent(firstDim);
  double asum = propensities(0);

  #ifdef OPTIMIZE_COMPUTE_PROPENSITIES
  while (asum <= atot2*r2)
  #else
  while (asum <= atot2*r2 && channel < nChannels-1)
    /* The condition channel < nChannels-1 is not necessary but is there to
       prevent numerical problems with the sum of propensities. For example,
       if the sum atot2 is not slightly inferior to the real sum, then we could
       have a channel number greater than the real number of channels. */
  #endif //OPTIMIZE_COMPUTE_PROPENSITIES
  {
    channel++;
    asum += propensities(channel);
  }

  ///- Get the time of the following division event in the whole cell collection.
  ///  The next division event will happen at time 't1'.
  double timeNextDivision = simulatorPtr_->cellCollection_.getTimeNextDivision();

  /**- If the time of next reaction is smaller than the time of next division
   *   event 't1', then
   *   - Apply reaction 'mu' in cell 'k' and update simulation time.
   * - If not, update simulation time to the time of next division
   *   event 't1' and divide the cell(s).
   */
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



