/***************************************************************************//**
 * Project: Colony
 *
 * \file    TimeSliceCellsAvg.h
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

#ifndef TIMESLICECELLSAVG_H
#define TIMESLICECELLSAVG_H

#include "compilation_options.h"

// standard C++ header files
#include <iostream>

// libraries header files
#include <blitz/array.h>

// user header files

// namespaces
using blitz::Array;
using std::cout;
using std::ostream;
using std::endl;



/**
 * Average of output values for one time slice.
 */
class TimeSliceCellsAvg
{
public:

  TimeSliceCellsAvg();
  ~TimeSliceCellsAvg();
  TimeSliceCellsAvg(const TimeSliceCellsAvg& t);

  void initialize(double timeStart, int stateAvgSize);

  void print ();

  TimeSliceCellsAvg& operator + (const TimeSliceCellsAvg& t);

  double timeStart_;
  Array<double,1> stateAvg_;
  double nCellsAvg_;
  int nSamples_;

};

#endif // TIMESLICECELLSAVG_H
