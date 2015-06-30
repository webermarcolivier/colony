/***************************************************************************//**
 * Project: Colony
 *
 * \file    TimeSlice.h
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

#ifndef TIMESLICE_H
#define TIMESLICE_H

#include "debug.h"

// standard C++ header files

// libraries header files
#include <blitz/array.h>

// user header files

// namespaces
using blitz::Array;
using blitz::TinyVector;
using blitz::firstDim;
using std::cout;
using std::endl;




/**
 * %Output values for one time slice.
 */
class TimeSlice
{
public:

  TimeSlice();
  ~TimeSlice();
  TimeSlice(const TimeSlice& t);

  void initialize
  (
    int nCells,
    double time,
    Array<double,1> volumeArray,
    Array<int,1>& x0Copy,
    Array<double,1> xConc,
    Array<TinyVector<double,3>,1> positionArray,
    Array<double,1> angleArray
  );

  void setPosition( Array<TinyVector<double,3>,1> &newPosition );
  void setAngle( Array<double,1> &newAngle );

  int nCells_;
  double time_;
  Array<double,1> volumeArray_;
  Array<double,1> stateAvg_;
  Array<int,1> t0State_;
  Array<double,1> xConc_;
  Array<TinyVector<double,3>,1> positionArray_;
  Array<double,1> angleArray_;
};

#endif // TIMESLICE_H
