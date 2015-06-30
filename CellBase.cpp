/***************************************************************************//**
 * Project: Colony
 *
 * \file    CellBase.cpp
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

#include "CellBase.h"

//------------------------------------------------------------------------------

CellBase::CellBase ()
{}

//------------------------------------------------------------------------------

CellBase::~CellBase ()
{}

//------------------------------------------------------------------------------

void CellBase::initialize (CellBaseInitParam& p)
{
  // Test if array 'x' has the right size.
  int nSpecies = p.chemicalSystemInitParam_.nSpecies_;
  int xSize = p.stateInitParam_.x_.size();
  if ( nSpecies != xSize )
  {
    cout << "ERROR: class CellBase, initialize(), x array does not have"
      " the correct size:\n";
    cout << "nSpecies = " << nSpecies << "\n";
    cout << "size(x0) = " << xSize << endl;
    exit(1);
  }

  // Test if 'stoichMatrix' has the right dimensions.
  int nrows = p.chemicalSystemInitParam_.stoichMatrix_.rows();
  int ncols = p.chemicalSystemInitParam_.stoichMatrix_.cols();
  int nChannels = p.chemicalSystemInitParam_.nChannels_;
  if (nrows != nChannels || ncols != nSpecies)
  {
    cout << "ERROR: class CellBase initialize(), stoichMatrix array does not have"
      " the correct size:\n";
    cout << "nChannels x nSpecies = " << nChannels << " x " << nSpecies << "\n";
    cout << "dim(stoichMatrix) = " << nrows << " x " << ncols << endl;
    exit(1);
  }

  // Set the size of propensities array in class State.
  p.stateInitParam_.propensitiesSize_ =
    p.chemicalSystemInitParam_.nChannels_;

  // Initialize State member.
  state_.initialize(p.stateInitParam_);

  // Initialize ChemicalSystem member.
  chemicalSystem_.initialize(p.chemicalSystemInitParam_);

}

//------------------------------------------------------------------------------

ostream& operator << (ostream &out, const CellBase& c)
{
  out << c.state_;
      //<< "nSpecies = " << c.chemicalSystem_.nSpecies_ << '\n'
      //<< "nChannels = " << c.chemicalSystem_.nChannels_ << '\n'
      /*<< "Dimensions(stoichMatrix) = "
        << c.chemicalSystem_.stoichMatrix_.extent(firstDim) << " x "
        << c.chemicalSystem_.stoichMatrix_.extent(secondDim) << '\n'*/
  out << endl;
  return out;
}

//------------------------------------------------------------------------------

void CellBase::changeVolume0By(double volume0Change)
{
  state_.cellVolume0_ += volume0Change;
}

//------------------------------------------------------------------------------

void CellBase::applyReaction(int mu)
{
  // Add the stoichiometric vector to the state.
  int i;
  for (i=0; i<chemicalSystem_.nSpecies_; i++)
  {
    state_.x_(i) += chemicalSystem_.stoichMatrix_(mu, i);
  }

  #ifdef X_POSITIVITY_CHECK
  checkPositivity();
  #endif //X_POSITIVITY_CHECK
}

//------------------------------------------------------------------------------

void CellBase::applyInterCellularReaction(Array<int,1> stoichVector)
{
  state_.x_ += stoichVector;

  #ifdef X_POSITIVITY_CHECK
  checkPositivity();
  #endif //X_POSITIVITY_CHECK
}

//------------------------------------------------------------------------------

void CellBase::computePropensities()
{
  chemicalSystem_.computePropensitiesPtr_(state_.x_,state_.propensities_);
}

//------------------------------------------------------------------------------

void CellBase::computePropensities(const double time)
{
  chemicalSystem_.computePropensitiesTimeDependentPtr_
  (
    state_.x_,
    getTimeDependentVolume(time),
    state_.cellVolume0_,
    state_.propensities_
  );
}

//------------------------------------------------------------------------------

void CellBase::computePropensitiesTimeDependentReactions(const double time)
{
  // To implement: a function written by Mathematica that compute only the
  // time-dependent reaction channels, for example only the 2nd order reactions.
  exit(1);
}

//------------------------------------------------------------------------------

void CellBase::setReferenceX(Array<int,1>& refArray)
{
  Array<int,1> arrayCopy( state_.x_.copy() );
  state_.x_.reference(refArray);
  state_.x_ = arrayCopy;
}

//------------------------------------------------------------------------------

void CellBase::setReferenceXConc(Array<double,1>& refArray)
{
  Array<double,1> arrayCopy( state_.xConc_.copy() );
  state_.xConc_.reference(refArray);
  state_.xConc_ = arrayCopy;
}

//------------------------------------------------------------------------------

void CellBase::setReferencePropensities(Array<double,1> refArray)
{
  Array<double,1> arrayCopy( state_.propensities_.copy() );
  state_.propensities_.reference(refArray);
  state_.propensities_ = arrayCopy;
}

//------------------------------------------------------------------------------

#ifdef X_POSITIVITY_CHECK
void CellBase::checkPositivity() const
{
  int i;
  for (i=0; i<(state_.x_).size(); i++)
  {
    if ( state_.x_(i) < 0 )
    {
      cout << "WARNING: class CellBase, function checkPositivity(),"
              " state vector x_ contains a negative number of molecules: "
              "x(" << i << ") = " << state_.x_(i) << endl;
    }
  }
}
#endif //X_POSITIVITY_CHECK

//------------------------------------------------------------------------------

void CellBase::setAngle(double angle)
{
  state_.angle_ = angle;
}

//------------------------------------------------------------------------------

void CellBase::setPosition(TinyVector<double,3> position)
{
  state_.position_ = position;
}

//------------------------------------------------------------------------------

void CellBase::copy(const CellBase& c2)
{
  state_ = c2.state_; // assigment operator of class State.
  chemicalSystem_ = c2.chemicalSystem_; // assigment operator of class ChemicalSystem.
}

//------------------------------------------------------------------------------

vector<string> CellBase::getListSpeciesName() const
{
  vector<string> listSpeciesName = state_.listSpeciesName_;
  return listSpeciesName;
}

//------------------------------------------------------------------------------

Array<double,1> CellBase::getXconc(const double time) const
{
  Array<double,1> xConc;
  Array<int,1> x ( state_.x_.copy() );
  xConc.resize( x.size() );
  double volume = getTimeDependentVolume(time);
  double conversionCoeffNM = 1.66054;
  int i;
  for(i=0; i<xConc.size(); ++i)
  {
    xConc(i) = conversionCoeffNM * double(x(i)) / volume;
  }
  return xConc;
}

//------------------------------------------------------------------------------

TinyVector<double,3> CellBase::getTimeDependentCellDimensions(const double time)
{
  return state_.getTimeDependentCellDimensions(time);
}

//------------------------------------------------------------------------------
