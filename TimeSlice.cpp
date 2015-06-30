/***************************************************************************//**
 * Project: Colony
 *
 * \file    TimeSlice.cpp
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

#include "TimeSlice.h"

//------------------------------------------------------------------------------

TimeSlice::TimeSlice()
{}

//------------------------------------------------------------------------------

TimeSlice::~TimeSlice()
{}

//------------------------------------------------------------------------------

TimeSlice::TimeSlice(const TimeSlice& t)
  : nCells_(t.nCells_),
    time_(t.time_),
    volumeArray_(t.volumeArray_.copy()),
    stateAvg_(t.stateAvg_.copy()),
    t0State_(t.t0State_.copy()),
    xConc_(t.xConc_.copy()),
    positionArray_(t.positionArray_.copy()),
    angleArray_(t.angleArray_.copy())
{}

//------------------------------------------------------------------------------

void TimeSlice::initialize
(
  int nCells,
  double time,
  Array<double,1> volumeArray,
  Array<int,1>& x0Copy,
  Array<double,1> xConc,
  Array<TinyVector<double,3>,1> positionArray,
  Array<double,1> angleArray
)
{
  nCells_ = nCells;
  time_   = time;
  volumeArray_.resize( volumeArray.size() );
  volumeArray_ = volumeArray;
  stateAvg_.resize( x0Copy.size() );
  stateAvg_ = 0.0;
  t0State_.resize( x0Copy.size() );
  t0State_ = x0Copy;
  xConc_.resize( xConc.size() );
  xConc_ = xConc;
  positionArray_.resize( positionArray.size() );
  positionArray_ = positionArray;
  angleArray_.resize( angleArray.size() );
  angleArray_ = angleArray;
}

//------------------------------------------------------------------------------

void TimeSlice::setPosition( Array<TinyVector<double,3>,1> &newPosition )
{
  positionArray_.resize( newPosition.size() );
  positionArray_ = newPosition;
}

//------------------------------------------------------------------------------

void TimeSlice::setAngle( Array<double,1> &newAngle )
{
  angleArray_.resize( newAngle.size() );
  angleArray_ = newAngle;
}

//------------------------------------------------------------------------------
