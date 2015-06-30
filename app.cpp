/***************************************************************************//**
 * Project: Colony
 *
 * \file    app.cpp
 * \author  Marc Weber\n
 *          The Si.M.Bio.Sys. Group (CosmoLab)\n
 *          Parc Científic de Barcelona\n
 *          Barcelona, Spain.\n
 *          http://thesimbiosys.isgreat.org
 * \version 0.1
 * \date    11/2009
 *
 *          Copyright 2009 by Marc Weber
 ******************************************************************************/

#include "debug.h"

// standard C++ header files
#include <iostream>
#include <iomanip>
#include <utility>

// libraries header files
#include <blitz/array.h>

// user header files
#include "Simulator.h"

// namespaces
using blitz::Array;
using blitz::Range;
using std::cout;
using std::ostream;
using std::endl;
using std::setw;



int main()
{

  Simulator simulator;

  //simulator.computeNTrajectories();

  simulator.computeParameterSet();

  return 0;
}
