/***************************************************************************//**
 * Project: Colony
 *
 * \file    TimeSliceCellsAvg.cpp
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

#include "TimeSliceCellsAvg.h"

//------------------------------------------------------------------------------

TimeSliceCellsAvg::TimeSliceCellsAvg()
  : nSamples_(0),
    timeStart_(-999.9)
{}

//------------------------------------------------------------------------------

TimeSliceCellsAvg::~TimeSliceCellsAvg()
{}

//------------------------------------------------------------------------------

TimeSliceCellsAvg::TimeSliceCellsAvg(const TimeSliceCellsAvg& t)
  : timeStart_(t.timeStart_),
    nSamples_(t.nSamples_),
    nCellsAvg_(t.nCellsAvg_),
    stateAvg_(t.stateAvg_.copy())
{}

//------------------------------------------------------------------------------

void TimeSliceCellsAvg::initialize(double timeStart, int stateAvgSize)
{
  timeStart_ = timeStart;
  stateAvg_.resize( stateAvgSize );
  stateAvg_ = 0.0;
  nSamples_ = 0;
  nCellsAvg_ = 0.0;
}

//------------------------------------------------------------------------------

TimeSliceCellsAvg& TimeSliceCellsAvg::operator + (const TimeSliceCellsAvg& t)
{
  // Test if the timeStart is the same in both objects.
  if ( timeStart_ != t.timeStart_ )
  {
    cout << "ERROR: class TimeSliceCellsAvg, operator +,"
            " the timeStart_ value is not the same in both slices." << endl;
    exit(1);
  }

  // Test if the stateAvg_ array have the same size.
  if ( stateAvg_.size() != t.stateAvg_.size() )
  {
    cout << "ERROR: class TimeSliceCellsAvg, operator +,"
            " arrays stateAvg_ do not have the same size.\n"
         << "size(stateAvg_) = " << stateAvg_.size()
         << " size(t.stateAvg_) = " << t.stateAvg_.size()
         << endl;
    exit(1);
  }

  stateAvg_ += t.stateAvg_;
  nCellsAvg_ += t.nCellsAvg_;
  nSamples_ ++;

  return *this;
}

//------------------------------------------------------------------------------

void TimeSliceCellsAvg::print ()
{

  Array<double,1> arrayCopy(stateAvg_.copy());
  //cout << "not normalized:\n" << arrayCopy << endl;
  arrayCopy /= double(nSamples_);
  cout << "normalized:\n" << arrayCopy << endl;
}

//------------------------------------------------------------------------------















