/***************************************************************************//**
 * Project: Colony
 *
 * \file    Input.h
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

#ifndef INPUT_H
#define INPUT_H

#include "debug.h"

// standard C++ header files
#include <iostream>   // input/output interface
#include <iomanip>    // input/output formatting
#include <sstream>    // manipulate strings as stream
#include <vector>
#include <string>     // manipulate strings of characters

// libraries header files
#include <blitz/array.h>

// user header files
#include "iotools.h"
#include "RandomNumberGenerator.h"
#include "Param/CellCollectionParam.h"
#include "Param/OutputParam.h"
#include "Param/SimulatorParam.h"

/**
 * Include the file for definition of the computePropensities functions.
 */
#include "computePropensitiesFunctions.h"

// namespaces
using blitz::Array;
using blitz::Range;
using blitz::firstDim;
using blitz::secondDim;
using blitz::thirdDim;
using std::cout;
using std::ostream;
using std::endl;
using std::setw;
using std::string;
using std::stringstream;
using std::vector;



/**
 * %Input parameters read from files "input.dat" and "input_simulation.dat"
 * in executable's directory.
 */
class Input
{
public:

  /**
   * The default constructor reads the input from file "input.dat" and
   * "input_simulation.dat"
   */
  Input(string inputFile);

  ~Input();

//---- OPERATIONS

  void setGlobalParameter(const double time);

  void setIndexParameterSet(const int i);

  int getNParameterSet() const;

  void setX0(const int iParameterSet);

//---- DATA

  int seed;

  static double ODEEquilibrationTime;

  CellCollectionParam p;

  SimulatorParam sim;

  OutputParam out;

  static Array<double,1> globalParameter;

  static Array<double,3> globalParameterTimePoints;

  static vector<string> listGlobalParametersName;


private:

  Array<double,1> simTimePoints_;

  int indexParameterSet_;

  int nParameterSet_;

  bool x0CellChangingWithGlobalParameter_;

  bool x0MilieuChangingWithGlobalParameter_;

  Array<int,2> x0CellTable_;

  bool x0CellHeterogenous_;

  Array<int,3> x0CellHeterogenousTable_;

  Array<double,2> initialPhaseTable_;

  Array<int,2> x0MilieuTable_;
};


#endif // INPUT_H
