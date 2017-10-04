/***************************************************************************//**
 * Project: Colony
 *
 * \file    computePropensitiesFunctions.h
 * \author  Marc Weber\n
 *          The Si.M.Bio.Sys. Group (CosmoLab)\n
 *          Parc Cienti­fic de Barcelona\n
 *          Barcelona, Spain.\n
 *          http://thesimbiosys.isgreat.org
 * \version 0.1
 * \date    03/2009
 *
 *          Output file from Mathematica.\n
 *          Copyright 2009 by Marc Weber
 ******************************************************************************/

#ifndef COMPUTEPROPENSITIESFUNCTIONS_CPP
#define COMPUTEPROPENSITIESFUNCTIONS_CPP

// libraries header files
#include <blitz/array.h>

// standard C++ header files
#include <iostream>

// libraries header files
#include <blitz/array.h>

// user header files
#include "Input.h"
#include "RandomNumberGenerator.h"

// namespaces
using blitz::Array;
using std::cout;
using std::endl;


void computePropensitiesCell
(
  const Array<int,1>& x,
  Array<double,1>& a
);


void computePropensitiesMilieu
(
  const Array<int,1>& x,
  Array<double,1>& a
);


void computePropensitiesCellTimeDependent
(
  const Array<int,1>& x,
  const double volume,
  const double volume0,
  Array<double,1>& a
);


void computePropensitiesMilieuTimeDependent
(
  const Array<int,1>& x,
  const double volume,
  const double volume0,
  Array<double,1>& a
);


void computePropensitiesCellMilieu
(
  const Array<int,1>& x1,
  const Array<int,1>& x2,
  const double volume1,
  const double volume2,
  Array<double,1>& a
);


void computePropensitiesTimeDependentDiffusion
(
  const Array<int,1>& x1,
  const Array<int,1>& x2,
  const double volume,
  const double volume0,
  const double volumeExt,
  Array<double,1>& a
);



#endif // COMPUTEPROPENSITIESFUNCTIONS_CPP

