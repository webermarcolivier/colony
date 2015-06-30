/***************************************************************************//**
 * Project: Colony
 *
 * \file    ChemicalLangevinComputeIncrementProfiling.h
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

#ifndef ChemicalLangevinComputeIncrementProfiling_H
#define ChemicalLangevinComputeIncrementProfiling_H

extern "C"
{
   #include <cblas.h>
}

// standard C++ header files
#include <iostream>

// libraries header files
#include <blitz/array.h>
//#include <armadillo>
//#include <cblas.h>

// user header files
#include "Input.h"
#include "RandomNumberGenerator.h"

// namespaces
using blitz::Array;
using blitz::Range;
using std::cout;
using std::endl;
//using arma::mat;
//using arma::randu;

/**
  * Algorithm of Euler for the chemical Langevin equation.
  * This is the same class as the #ChemicalLangevinComputeIncrement but
  * I separated all the different parts in the algorithm into different
  * private methods in order to be able to profile the code.
  * The profiler will tell us how much time the algorithm
  * spends in each part (method) of the algorithm. Note: the line by line
  * profiling does not work with the Blitz++ library we use because
  * of debug options.
  */
class ChemicalLangevinComputeIncrementProfiling
{
public:

  ChemicalLangevinComputeIncrementProfiling();
  ~ChemicalLangevinComputeIncrementProfiling();

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

  double computeMean(Array<double,1>& array);

  void computeDeterministicIncrement
  (
    Array<double,1>& xConcCells,
    Array<double,1>& xConcMilieu,
    const int nCells,
    const double volumeCell,
    const double volumeExt,
    const double chemicalLangevinTimeStep,
    double avgU,
    double avgV
  );

  void computeBNoise
  (
    Array<double,2>& bNoise,
    Array<double,1>& xConcCells,
    Array<double,1>& xConcMilieu,
    double avgU,
    double avgV
  );

  void generateRandomVector
  (
    Array<double,1>& randomGaussianVector
  );

  void computeMatrixVectorProduct
  (
    Array<double,2>& bNoise,
    Array<double,1>& randomGaussianVector,
    Array<double,1>& noiseIncrementVector
  );

  void addNoiseIncrement
  (
    Array<double,1>& xConcCells,
    Array<double,1>& xConcMilieu,
    const int nCells,
    const double volumeCell,
    const double volumeExt,
    const double chemicalLangevinTimeStep
  );


};

#endif // ChemicalLangevinComputeIncrementProfiling_H
