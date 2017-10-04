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

#include "compilation_options.h"

#ifdef USE_CHEMICAL_LANGEVIN

#include "chemicalLangevinComputeIncrement.h"

//------------------------------------------------------------------------------

ChemicalLangevinComputeIncrement::ChemicalLangevinComputeIncrement()
{
  nSCell = 7;
  nSMilieu = 1;
  nRCell = 26;
  nRMilieu = 4;
  nRCellMilieu = 2;

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
  //double avgU = mean( xConcCells( blitz::Range( -1 + Extract[{}, {1, 1}], blitz::toEnd, nSCell ) ) );
  //double avgV = mean( xConcCells( blitz::Range( -1 + Extract[{}, {1, 1}], blitz::toEnd, nSCell ) ) );

 

 // Computing deterministic increments

  // Compute increments for cell species
  for (int iCell=0; iCell<nCells; iCell++)
  {
    xConcCells(iCell*nSCell + 0) += (0.1*xConcCells(5 + 7*iCell) + 10.*xConcCells(6 + 7*iCell) - 0.01586294361119891*xConcCells(7*iCell)) * chemicalLangevinTimeStep;
    xConcCells(iCell*nSCell + 1) += (0.01*xConcCells(5 + 7*iCell) + xConcCells(6 + 7*iCell) - 0.01586294361119891*xConcCells(1 + 7*iCell) - 0.1*xConcCells(2 + 7*iCell)*xConcCells(1 + 7*iCell) + 10.*xConcCells(3 + 7*iCell)) * chemicalLangevinTimeStep;
    xConcCells(iCell*nSCell + 2) += (-10.0158629436112*xConcCells(2 + 7*iCell) + 10.*xConcMilieu(0) + 0.04*xConcCells(7*iCell) - 0.1*xConcCells(2 + 7*iCell)*xConcCells(1 + 7*iCell) + 10.*xConcCells(3 + 7*iCell)) * chemicalLangevinTimeStep;
    xConcCells(iCell*nSCell + 3) += (0.1*xConcCells(2 + 7*iCell)*xConcCells(1 + 7*iCell) - 10.0158629436112*xConcCells(3 + 7*iCell) - 0.1*pow(xConcCells(3 + 7*iCell),2) + 2.*xConcCells(4 + 7*iCell)) * chemicalLangevinTimeStep;
    xConcCells(iCell*nSCell + 4) += (10.*xConcCells(6 + 7*iCell) + 0.05*pow(xConcCells(3 + 7*iCell),2) - 1.015862943611199*xConcCells(4 + 7*iCell) - 0.05*xConcCells(5 + 7*iCell)*xConcCells(4 + 7*iCell)) * chemicalLangevinTimeStep;
    xConcCells(iCell*nSCell + 5) += (10.0138629436112*xConcCells(6 + 7*iCell) - 0.05*xConcCells(5 + 7*iCell)*xConcCells(4 + 7*iCell)) * chemicalLangevinTimeStep;
    xConcCells(iCell*nSCell + 6) += (-10.0138629436112*xConcCells(6 + 7*iCell) + 0.05*xConcCells(5 + 7*iCell)*xConcCells(4 + 7*iCell)) * chemicalLangevinTimeStep;
  }

  // Compute increments for milieu species. We need to compute separately the increments
  // for the milieu reactions and for the cell-milieu reactions, which will be corrected for volume effects.

  // Compute increments for milieu species for milieu reactions.
  xConcMilieu(0) += (-0.002*xConcMilieu(0)) * chemicalLangevinTimeStep;

  // Compute increments for milieu species for cell-milieu reactions.
  for (int iCell=0; iCell<nCells; iCell++)
  {
    xConcMilieu(0) += ((10.*volumeCell*xConcCells(2 + 7*iCell))/volumeExt - (10.*volumeCell*xConcMilieu(0))/volumeExt) * chemicalLangevinTimeStep;
  }



  // Computing stochastic increments


  // Compute Bnoise matrix
  //double noiseTermDAvgU = sqrt(avgU*Dlattice);
  //double noiseTermDAvgV = sqrt(avgV*Dlattice);

  bNoise(0,0) = Bnoise[1][1];
  bNoise(0,1) = Bnoise[1][2];
  bNoise(0,2) = Bnoise[1][3];
  bNoise(0,3) = Bnoise[1][4];
  bNoise(0,4) = Bnoise[1][5];
  bNoise(0,5) = Bnoise[1][6];
  bNoise(0,6) = Bnoise[1][7];
  bNoise(0,7) = Bnoise[1][8];
  bNoise(0,8) = Bnoise[1][9];
  bNoise(0,9) = Bnoise[1][10];
  bNoise(0,10) = Bnoise[1][11];
  bNoise(0,11) = Bnoise[1][12];
  bNoise(0,12) = Bnoise[1][13];
  bNoise(0,13) = Bnoise[1][14];
  bNoise(0,14) = Bnoise[1][15];
  bNoise(0,15) = Bnoise[1][16];
  bNoise(0,16) = Bnoise[1][17];
  bNoise(0,17) = Bnoise[1][18];
  bNoise(0,18) = Bnoise[1][19];
  bNoise(0,19) = Bnoise[1][20];
  bNoise(0,20) = Bnoise[1][21];
  bNoise(0,21) = Bnoise[1][22];
  bNoise(0,22) = Bnoise[1][23];
  bNoise(0,23) = Bnoise[1][24];
  bNoise(0,24) = Bnoise[1][25];
  bNoise(0,25) = Bnoise[1][26];
  bNoise(0,26) = Bnoise[1][27];
  bNoise(0,27) = Bnoise[1][28];
  bNoise(0,28) = Bnoise[1][29];
  bNoise(0,29) = Bnoise[1][30];
  bNoise(0,30) = Bnoise[1][31];
  bNoise(0,31) = Bnoise[1][32];
  bNoise(1,0) = Bnoise[2][1];
  bNoise(1,1) = Bnoise[2][2];
  bNoise(1,2) = Bnoise[2][3];
  bNoise(1,3) = Bnoise[2][4];
  bNoise(1,4) = Bnoise[2][5];
  bNoise(1,5) = Bnoise[2][6];
  bNoise(1,6) = Bnoise[2][7];
  bNoise(1,7) = Bnoise[2][8];
  bNoise(1,8) = Bnoise[2][9];
  bNoise(1,9) = Bnoise[2][10];
  bNoise(1,10) = Bnoise[2][11];
  bNoise(1,11) = Bnoise[2][12];
  bNoise(1,12) = Bnoise[2][13];
  bNoise(1,13) = Bnoise[2][14];
  bNoise(1,14) = Bnoise[2][15];
  bNoise(1,15) = Bnoise[2][16];
  bNoise(1,16) = Bnoise[2][17];
  bNoise(1,17) = Bnoise[2][18];
  bNoise(1,18) = Bnoise[2][19];
  bNoise(1,19) = Bnoise[2][20];
  bNoise(1,20) = Bnoise[2][21];
  bNoise(1,21) = Bnoise[2][22];
  bNoise(1,22) = Bnoise[2][23];
  bNoise(1,23) = Bnoise[2][24];
  bNoise(1,24) = Bnoise[2][25];
  bNoise(1,25) = Bnoise[2][26];
  bNoise(1,26) = Bnoise[2][27];
  bNoise(1,27) = Bnoise[2][28];
  bNoise(1,28) = Bnoise[2][29];
  bNoise(1,29) = Bnoise[2][30];
  bNoise(1,30) = Bnoise[2][31];
  bNoise(1,31) = Bnoise[2][32];
  bNoise(2,0) = Bnoise[3][1];
  bNoise(2,1) = Bnoise[3][2];
  bNoise(2,2) = Bnoise[3][3];
  bNoise(2,3) = Bnoise[3][4];
  bNoise(2,4) = Bnoise[3][5];
  bNoise(2,5) = Bnoise[3][6];
  bNoise(2,6) = Bnoise[3][7];
  bNoise(2,7) = Bnoise[3][8];
  bNoise(2,8) = Bnoise[3][9];
  bNoise(2,9) = Bnoise[3][10];
  bNoise(2,10) = Bnoise[3][11];
  bNoise(2,11) = Bnoise[3][12];
  bNoise(2,12) = Bnoise[3][13];
  bNoise(2,13) = Bnoise[3][14];
  bNoise(2,14) = Bnoise[3][15];
  bNoise(2,15) = Bnoise[3][16];
  bNoise(2,16) = Bnoise[3][17];
  bNoise(2,17) = Bnoise[3][18];
  bNoise(2,18) = Bnoise[3][19];
  bNoise(2,19) = Bnoise[3][20];
  bNoise(2,20) = Bnoise[3][21];
  bNoise(2,21) = Bnoise[3][22];
  bNoise(2,22) = Bnoise[3][23];
  bNoise(2,23) = Bnoise[3][24];
  bNoise(2,24) = Bnoise[3][25];
  bNoise(2,25) = Bnoise[3][26];
  bNoise(2,26) = Bnoise[3][27];
  bNoise(2,27) = Bnoise[3][28];
  bNoise(2,28) = Bnoise[3][29];
  bNoise(2,29) = Bnoise[3][30];
  bNoise(2,30) = Bnoise[3][31];
  bNoise(2,31) = Bnoise[3][32];
  bNoise(3,0) = Bnoise[4][1];
  bNoise(3,1) = Bnoise[4][2];
  bNoise(3,2) = Bnoise[4][3];
  bNoise(3,3) = Bnoise[4][4];
  bNoise(3,4) = Bnoise[4][5];
  bNoise(3,5) = Bnoise[4][6];
  bNoise(3,6) = Bnoise[4][7];
  bNoise(3,7) = Bnoise[4][8];
  bNoise(3,8) = Bnoise[4][9];
  bNoise(3,9) = Bnoise[4][10];
  bNoise(3,10) = Bnoise[4][11];
  bNoise(3,11) = Bnoise[4][12];
  bNoise(3,12) = Bnoise[4][13];
  bNoise(3,13) = Bnoise[4][14];
  bNoise(3,14) = Bnoise[4][15];
  bNoise(3,15) = Bnoise[4][16];
  bNoise(3,16) = Bnoise[4][17];
  bNoise(3,17) = Bnoise[4][18];
  bNoise(3,18) = Bnoise[4][19];
  bNoise(3,19) = Bnoise[4][20];
  bNoise(3,20) = Bnoise[4][21];
  bNoise(3,21) = Bnoise[4][22];
  bNoise(3,22) = Bnoise[4][23];
  bNoise(3,23) = Bnoise[4][24];
  bNoise(3,24) = Bnoise[4][25];
  bNoise(3,25) = Bnoise[4][26];
  bNoise(3,26) = Bnoise[4][27];
  bNoise(3,27) = Bnoise[4][28];
  bNoise(3,28) = Bnoise[4][29];
  bNoise(3,29) = Bnoise[4][30];
  bNoise(3,30) = Bnoise[4][31];
  bNoise(3,31) = Bnoise[4][32];
  bNoise(4,0) = Bnoise[5][1];
  bNoise(4,1) = Bnoise[5][2];
  bNoise(4,2) = Bnoise[5][3];
  bNoise(4,3) = Bnoise[5][4];
  bNoise(4,4) = Bnoise[5][5];
  bNoise(4,5) = Bnoise[5][6];
  bNoise(4,6) = Bnoise[5][7];
  bNoise(4,7) = Bnoise[5][8];
  bNoise(4,8) = Bnoise[5][9];
  bNoise(4,9) = Bnoise[5][10];
  bNoise(4,10) = Bnoise[5][11];
  bNoise(4,11) = Bnoise[5][12];
  bNoise(4,12) = Bnoise[5][13];
  bNoise(4,13) = Bnoise[5][14];
  bNoise(4,14) = Bnoise[5][15];
  bNoise(4,15) = Bnoise[5][16];
  bNoise(4,16) = Bnoise[5][17];
  bNoise(4,17) = Bnoise[5][18];
  bNoise(4,18) = Bnoise[5][19];
  bNoise(4,19) = Bnoise[5][20];
  bNoise(4,20) = Bnoise[5][21];
  bNoise(4,21) = Bnoise[5][22];
  bNoise(4,22) = Bnoise[5][23];
  bNoise(4,23) = Bnoise[5][24];
  bNoise(4,24) = Bnoise[5][25];
  bNoise(4,25) = Bnoise[5][26];
  bNoise(4,26) = Bnoise[5][27];
  bNoise(4,27) = Bnoise[5][28];
  bNoise(4,28) = Bnoise[5][29];
  bNoise(4,29) = Bnoise[5][30];
  bNoise(4,30) = Bnoise[5][31];
  bNoise(4,31) = Bnoise[5][32];
  bNoise(5,0) = Bnoise[6][1];
  bNoise(5,1) = Bnoise[6][2];
  bNoise(5,2) = Bnoise[6][3];
  bNoise(5,3) = Bnoise[6][4];
  bNoise(5,4) = Bnoise[6][5];
  bNoise(5,5) = Bnoise[6][6];
  bNoise(5,6) = Bnoise[6][7];
  bNoise(5,7) = Bnoise[6][8];
  bNoise(5,8) = Bnoise[6][9];
  bNoise(5,9) = Bnoise[6][10];
  bNoise(5,10) = Bnoise[6][11];
  bNoise(5,11) = Bnoise[6][12];
  bNoise(5,12) = Bnoise[6][13];
  bNoise(5,13) = Bnoise[6][14];
  bNoise(5,14) = Bnoise[6][15];
  bNoise(5,15) = Bnoise[6][16];
  bNoise(5,16) = Bnoise[6][17];
  bNoise(5,17) = Bnoise[6][18];
  bNoise(5,18) = Bnoise[6][19];
  bNoise(5,19) = Bnoise[6][20];
  bNoise(5,20) = Bnoise[6][21];
  bNoise(5,21) = Bnoise[6][22];
  bNoise(5,22) = Bnoise[6][23];
  bNoise(5,23) = Bnoise[6][24];
  bNoise(5,24) = Bnoise[6][25];
  bNoise(5,25) = Bnoise[6][26];
  bNoise(5,26) = Bnoise[6][27];
  bNoise(5,27) = Bnoise[6][28];
  bNoise(5,28) = Bnoise[6][29];
  bNoise(5,29) = Bnoise[6][30];
  bNoise(5,30) = Bnoise[6][31];
  bNoise(5,31) = Bnoise[6][32];
  bNoise(6,0) = Bnoise[7][1];
  bNoise(6,1) = Bnoise[7][2];
  bNoise(6,2) = Bnoise[7][3];
  bNoise(6,3) = Bnoise[7][4];
  bNoise(6,4) = Bnoise[7][5];
  bNoise(6,5) = Bnoise[7][6];
  bNoise(6,6) = Bnoise[7][7];
  bNoise(6,7) = Bnoise[7][8];
  bNoise(6,8) = Bnoise[7][9];
  bNoise(6,9) = Bnoise[7][10];
  bNoise(6,10) = Bnoise[7][11];
  bNoise(6,11) = Bnoise[7][12];
  bNoise(6,12) = Bnoise[7][13];
  bNoise(6,13) = Bnoise[7][14];
  bNoise(6,14) = Bnoise[7][15];
  bNoise(6,15) = Bnoise[7][16];
  bNoise(6,16) = Bnoise[7][17];
  bNoise(6,17) = Bnoise[7][18];
  bNoise(6,18) = Bnoise[7][19];
  bNoise(6,19) = Bnoise[7][20];
  bNoise(6,20) = Bnoise[7][21];
  bNoise(6,21) = Bnoise[7][22];
  bNoise(6,22) = Bnoise[7][23];
  bNoise(6,23) = Bnoise[7][24];
  bNoise(6,24) = Bnoise[7][25];
  bNoise(6,25) = Bnoise[7][26];
  bNoise(6,26) = Bnoise[7][27];
  bNoise(6,27) = Bnoise[7][28];
  bNoise(6,28) = Bnoise[7][29];
  bNoise(6,29) = Bnoise[7][30];
  bNoise(6,30) = Bnoise[7][31];
  bNoise(6,31) = Bnoise[7][32];
  bNoise(7,0) = Bnoise[8][1];
  bNoise(7,1) = Bnoise[8][2];
  bNoise(7,2) = Bnoise[8][3];
  bNoise(7,3) = Bnoise[8][4];
  bNoise(7,4) = Bnoise[8][5];
  bNoise(7,5) = Bnoise[8][6];
  bNoise(7,6) = Bnoise[8][7];
  bNoise(7,7) = Bnoise[8][8];
  bNoise(7,8) = Bnoise[8][9];
  bNoise(7,9) = Bnoise[8][10];
  bNoise(7,10) = Bnoise[8][11];
  bNoise(7,11) = Bnoise[8][12];
  bNoise(7,12) = Bnoise[8][13];
  bNoise(7,13) = Bnoise[8][14];
  bNoise(7,14) = Bnoise[8][15];
  bNoise(7,15) = Bnoise[8][16];
  bNoise(7,16) = Bnoise[8][17];
  bNoise(7,17) = Bnoise[8][18];
  bNoise(7,18) = Bnoise[8][19];
  bNoise(7,19) = Bnoise[8][20];
  bNoise(7,20) = Bnoise[8][21];
  bNoise(7,21) = Bnoise[8][22];
  bNoise(7,22) = Bnoise[8][23];
  bNoise(7,23) = Bnoise[8][24];
  bNoise(7,24) = Bnoise[8][25];
  bNoise(7,25) = Bnoise[8][26];
  bNoise(7,26) = Bnoise[8][27];
  bNoise(7,27) = Bnoise[8][28];
  bNoise(7,28) = Bnoise[8][29];
  bNoise(7,29) = Bnoise[8][30];
  bNoise(7,30) = Bnoise[8][31];
  bNoise(7,31) = Bnoise[8][32];


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
