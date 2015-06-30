/***************************************************************************//**
 * Project: Colony
 *
 * \file    State.cpp
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

#include "State.h"

//------------------------------------------------------------------------------

State::State()
  : x_(),  //Rem: the Array is 0-initialized, its data cannot be accessed
    xConc_(),
    propensities_(), //Rem: the Array is 0-initialized, its data cannot be accessed
    position_(0.0,0.0,0.0),
    angle_(0.0),
    cellLength_(0),
    cellLength0_(0),
    cellHeight_(0),
    cellVolume0_(0)
{
  timePreviousDivision_ = - numeric_limits<double>::infinity();
  cellCycleDuration_    =   numeric_limits<double>::infinity();
}

//------------------------------------------------------------------------------

State::~State()
{
  x_.free();
  xConc_.free();
  propensities_.free();
}

//------------------------------------------------------------------------------

void State::initialize(const StateInitParam& p)
{
  int size = p.x_.size();
  if (size > 0)
  {
    x_.resize(size);
    x_ = p.x_;
    xConc_.resize(size);
    double conversionCoeffNM = 1.66054;
    for(int i=0; i<xConc_.size(); ++i)
    {
      xConc_(i) = conversionCoeffNM * double(x_(i)) / p.cellVolume0_;
    }
  }

  listSpeciesName_ = p.listSpeciesName_;

  propensities_.resize( p.propensitiesSize_ );
  propensities_ = 0.0;

  cellVolume0_ = p.cellVolume0_;
  cellLength0_ = p.cellLength0_;
  cellHeight0_ = p.cellHeight0_;
  cellHeight_ = cellHeight0_;
  cellLength_ = cellLength0_;

  position_(0) = 0.0;
  position_(1) = 0.0;
  position_(2) = 0.0;
  angle_ = 0.0;
}

//------------------------------------------------------------------------------

void State::copy(const State& s2)
{
  int size = s2.x_.size();
  if (size > 0)
  {
    x_.resize(size);
    x_ = s2.x_;
  }

  size = s2.xConc_.size();
  if (size > 0)
  {
    xConc_.resize(size);
    xConc_ = s2.xConc_;
  }

  listSpeciesName_ = s2.listSpeciesName_;

  propensities_.resize( s2.propensities_.size() );
  propensities_ = s2.propensities_;

  cellVolume0_ = s2.cellVolume0_;
  cellLength0_ = s2.cellLength0_;
  cellHeight0_ = s2.cellHeight0_;
  cellHeight_ = s2.cellHeight_;
  cellLength_ = s2.cellLength_;

  timePreviousDivision_ = s2.timePreviousDivision_;
  cellCycleDuration_    = s2.cellCycleDuration_;

  position_ = s2.position_;
  angle_ = s2.angle_;
}

//------------------------------------------------------------------------------

State& State::operator = (const State& s2)
{
  // Test for self-assignment
  if (&s2 == this)
  {
    // Do nothing and just return myself.
    return *this;
  }

  copy(s2);

  return *this;
}

//------------------------------------------------------------------------------

ostream& operator<< (ostream& out, const State& s)
{
  out << "cellVolume0 = " << s.cellVolume0_ << '\n'
      << "position = " << "(" << s.position_(0) << "," << s.position_(1) << ")" << '\n'
      << "angle = " << s.angle_ << '\n'
      << "timePreviousDivision = " << s.timePreviousDivision_ << '\n'
      << "cellCycleDuration = " << s.cellCycleDuration_ << '\n'
      << "timeNextDivision = " << s.timePreviousDivision_ +
          s.cellCycleDuration_ << '\n';
  out << "x: " << '\n';
  int i;
  int minWidth = 4;
  vector<int> widthTable;
  for (i=0; i<s.x_.size(); i++)
  {
    out << setw(minWidth) << s.listSpeciesName_[i] << " ";
    int width = max( minWidth, (int)s.listSpeciesName_[i].length() );
    widthTable.push_back(width);
  }
  out << '\n';
  for (i=0; i<s.x_.size(); i++)
  {
    out << setw(widthTable[i]) << s.x_(i) << " ";
  }
  out << endl;
  return out;
}

//------------------------------------------------------------------------------

double State::getTimeDependentVolume(const double time) const
{
  #ifdef TIME_DEPENDENT_PROPENSITIES
    return cellVolume0_*pow(2.0,(time-timePreviousDivision_)/ cellCycleDuration_);
  #else
    return cellVolume0_;
  #endif
}

//------------------------------------------------------------------------------

TinyVector<double,3> State::getTimeDependentCellDimensions(const double time)
{
  TinyVector<double,3> dimensions;
  dimensions(0) = getTimeDependentVolume(time);
  dimensions(1) = cellHeight_;
  cellLength_ = ( getTimeDependentVolume(time) / cellVolume0_ ) * cellLength0_;
  dimensions(2) = cellLength_;

  return dimensions;
}

//------------------------------------------------------------------------------



