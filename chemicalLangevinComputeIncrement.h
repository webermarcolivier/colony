/***************************************************************************//**
 * Project: Colony
 *
 * \file    chemicalLangevinComputeIncrement.h
 * \author  Marc Weber\n
 *          The Si.M.Bio.Sys. Group (CosmoLab)\n
 *          Parc Cienti­fic de Barcelona\n
 *          Barcelona, Spain.\n
 * \version 0.1
 * \date    04/2012
 *
 *          Output file from Mathematica.\n
 *          Copyright 2009 by Marc Weber
 ******************************************************************************/

#ifndef CHEMICALLANGEVINCOMPUTEINCREMENT_H
#define CHEMICALLANGEVINCOMPUTEINCREMENT_H

extern "C"
{
   #include <cblas.h>
}

// standard C++ header files
#include <iostream>

// libraries header files
#include <blitz/array.h>

// user header files
#include "Input.h"
#include "RandomNumberGenerator.h"

// namespaces
using blitz::Array;
using blitz::Range;
using std::cout;
using std::endl;

/**
  * Algorithm of Euler for the chemical Langevin equation.
  * This class is written by the mathematica notebook.
  */
class ChemicalLangevinComputeIncrement
{
public:

  ChemicalLangevinComputeIncrement();
  ~ChemicalLangevinComputeIncrement();

void chemicalLangevinComputeIncrement
(
  Array<double,1>& xConcCells,
  Array<double,1>& xConcMilieu,
  const int nCells,
  const double volumeCell,
  const double volumeExt,
  const double chemicalLangevinTimeStep
);


  void setFixedNCells(int nCells);
  double proportionNegativeSteps;


private:

  Array<double,2> bNoise;
  Array<double,1> noiseIncrementVector;
  int nSCell;
  int nSMilieu;
  int nRCell;
  int nRMilieu;
  int nRCellMilieu;
  int nSCellstot;
  int nRtot;
  int nStot;
  double nMConversion;
  double nMConversionPower3;
  double sqrtDx;
  double sqrtTimeStep;

  double totalNumberSteps_;
  double nbStepsWithNegativeCondition_;

};

#endif // CHEMICALLANGEVINCOMPUTEINCREMENT_H
