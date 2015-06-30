/***************************************************************************//**
 * Project: Colony
 *
 * \file    chemicalLangevinComputeIncrement.cpp
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

#include "debug.h"

#ifdef USE_CHEMICAL_LANGEVIN

#include "chemicalLangevinComputeIncrement.h"

//------------------------------------------------------------------------------

ChemicalLangevinComputeIncrement::ChemicalLangevinComputeIncrement()
{
  nSCell = 2;
  nSMilieu = 2;
  nRCell = 4;
  nRMilieu = 4;
  nRCellMilieu = 4;

  // Conversion coefficient between number of molecules and concentration in nM.
  nMConversion = 1.6605392767355127;
  nMConversionPower3 = pow(nMConversion,3);

}

//------------------------------------------------------------------------------

ChemicalLangevinComputeIncrement::~ChemicalLangevinComputeIncrement()
{}

//------------------------------------------------------------------------------

void ChemicalLangevinComputeIncrement::setFixedNCells(int nCells)
{
  nSCellstot = nCells*nSCell;
  nRtot = 2.*(nCells*(nRCell+nRCellMilieu) + nRMilieu);
  nStot = nCells*nSCell + nSMilieu;

  bNoise.resize(nStot,nRtot);
  bNoise = 0.;

  noiseIncrementVector.resize(nStot);
}

//------------------------------------------------------------------------------

void ChemicalLangevinComputeIncrement::chemicalLangevinComputeIncrement
(
  Array<double,1>& xConcCells,
  Array<double,1>& xConcMilieu,
  const int nCells,
  const double volumeCell,
  const double volumeExt,
  const double chemicalLangevinTimeStep
)
{

  int nSCellstot = nCells*nSCell;
  nRtot = 2.*(nCells*(nRCell+nRCellMilieu) + nRMilieu);
  nStot = nCells*nSCell + nSMilieu;
  double sqrtTimeStep = sqrt(chemicalLangevinTimeStep);

  // Computing the population average of U and V
  //double avgU = mean( xConcCells( blitz::Range( 0, blitz::toEnd, nSCell ) ) );
  //double avgV = mean( xConcCells( blitz::Range( 1, blitz::toEnd, nSCell ) ) );

 
 // Computing the expression of alpha and beta in function of the diffusion coefficient D
  double alpha = (0.5000000000000007 + Input::globalParameter(0)*(1.500000000000001 + 1.0000000000000013*Input::globalParameter(0)))/pow(1. + 1.*Input::globalParameter(0),2);
  double beta = 20.000000000000004 - 10.000000000000002/(1. + 1.*Input::globalParameter(0));

 

 // Computing deterministic increments

  // Compute increments for cell species
  for (int iCell=0; iCell<nCells; iCell++)
  {
    xConcCells(iCell*nSCell + 0) += (alpha - 1.*xConcCells(2*iCell) - 1.*Input::globalParameter(0)*xConcCells(2*iCell) + Input::globalParameter(0)*xConcMilieu(0) + (0.125*beta)/(0.125 + xConcCells(1 + 2*iCell)*xConcCells(1 + 2*iCell)*xConcCells(1 + 2*iCell))) * chemicalLangevinTimeStep;
    xConcCells(iCell*nSCell + 1) += (alpha + (0.125*beta)/(0.125 + xConcCells(2*iCell)*xConcCells(2*iCell)*xConcCells(2*iCell)) - 1.*xConcCells(1 + 2*iCell) - 1.*Input::globalParameter(0)*xConcCells(1 + 2*iCell) + Input::globalParameter(0)*xConcMilieu(1)) * chemicalLangevinTimeStep;
  }

  // Compute increments for milieu species. We need to compute separately the increments
  // for the milieu reactions and for the cell-milieu reactions, which will be corrected for volume effects.

  // Compute increments for milieu species for milieu reactions.
  xConcMilieu(0) += (-1.*xConcMilieu(0)) * chemicalLangevinTimeStep;
  xConcMilieu(1) += (-1.*xConcMilieu(1)) * chemicalLangevinTimeStep;

  // Compute increments for milieu species for cell-milieu reactions.
  for (int iCell=0; iCell<nCells; iCell++)
  {
    xConcMilieu(0) += ((volumeCell*Input::globalParameter(0)*xConcCells(2*iCell))/volumeExt - (1.*volumeCell*Input::globalParameter(0)*xConcMilieu(0))/volumeExt) * chemicalLangevinTimeStep;
    xConcMilieu(1) += ((volumeCell*Input::globalParameter(0)*xConcCells(1 + 2*iCell))/volumeExt - (1.*volumeCell*Input::globalParameter(0)*xConcMilieu(1))/volumeExt) * chemicalLangevinTimeStep;
  }



  // Computing stochastic increments


  // Compute Bnoise matrix
  //double noiseTermDAvgU = sqrt(avgU*Dlattice);
  //double noiseTermDAvgV = sqrt(avgV*Dlattice);

  bNoise(0,0) = List(0.,0.);
  bNoise(0,3) = List(0.,0.);
  bNoise(0,4) = List(0.,0.);
  bNoise(0,7) = List(0.,0.);
  bNoise(1,0) = List(0.,0.);
  bNoise(1,3) = List(0.,0.);
  bNoise(1,4) = List(0.,0.);
  bNoise(1,7) = List(0.,0.);
  bNoise(2,0) = List(0.,0.);
  bNoise(2,3) = List(0.,0.);
  bNoise(2,4) = List(0.,0.);
  bNoise(2,7) = List(0.,0.);
  bNoise(3,0) = List(0.,0.);
  bNoise(3,3) = List(0.,0.);
  bNoise(3,4) = List(0.,0.);
  bNoise(3,7) = List(0.,0.);


  // DoWhile loop with condition on positive values of xConc
  bool thereIsANegativeValue = false;
  do 
  {

    // Generating a vector of nRtot gaussian normal random numbers
    Array<double,1> randomGaussianVector(nRtot);
    for (int i=0; i<nRtot; i++)
    {
      randomGaussianVector(i) = RandomNumberGenerator::getNormal();
    }

    // Compute the matrix-vector product noiseIncrementVector = BNoise * randomGaussianVector.
    // BLAS matrix-vector multiplication
    // Note: We use the optimized BLAS library from ATLAS
    // (Automatically Tuned Linear Algebra Software) which is
    // architecture-specific.
    // See http://math-atlas.sourceforge.net/

    const int m = bNoise.extent(blitz::firstDim);
    const int n = bNoise.extent(blitz::secondDim);
    const double alphaBLAS = 1.0;
    const double betaBLAS = 0.0;
    const double *bNoisePointer = bNoise.data();
    double *noiseIncrementVectorPointer = noiseIncrementVector.data();
    const double *randomGaussianVectorPointer = randomGaussianVector.data();
    const int lda = n;
    const int incX = 1;
    const int incY = 1;

    cblas_dgemv(CblasRowMajor,CblasNoTrans,m,n,alphaBLAS,bNoisePointer,lda,
                randomGaussianVectorPointer,incX,betaBLAS,
                noiseIncrementVectorPointer,incY);


    // Compute if the stochastic increment would result in a negative value of xConc
    thereIsANegativeValue = any( 
      xConcCells + (sqrt(nMConversion/volumeCell))*noiseIncrementVector(Range(0,nSCellstot-1))*sqrtTimeStep < 0 )
      || 
      any( xConcMilieu  + (sqrt(nMConversion/volumeExt))*noiseIncrementVector(Range(nSCellstot,nStot-1))*sqrtTimeStep < 0 );

    if(thereIsANegativeValue) nbStepsWithNegativeCondition_ = nbStepsWithNegativeCondition_ + 1.0;
    totalNumberSteps_ = totalNumberSteps_ + 1.0;
  
  } while (thereIsANegativeValue);

  // Now that we know that the stochastic increment does not lead to a negative value, we apply the increment.
  xConcCells += (sqrt(nMConversion/volumeCell))*noiseIncrementVector(Range(0,nSCellstot-1))*sqrtTimeStep;
  xConcMilieu += (sqrt(nMConversion/volumeExt))*noiseIncrementVector(Range(nSCellstot,nStot-1))*sqrtTimeStep;

  proportionNegativeSteps = nbStepsWithNegativeCondition_ / totalNumberSteps_;

}

#endif // USE_CHEMICAL_LANGEVIN
