/***************************************************************************//**
 * Project: Colony
 *
 * \file    ChemicalLangevinComputeIncrementProfiling.cpp
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

#include "chemicalLangevinComputeIncrementProfiling.h"

//------------------------------------------------------------------------------

ChemicalLangevinComputeIncrementProfiling::ChemicalLangevinComputeIncrementProfiling()
{
  nSCell = 2;
  nSMilieu = 1;
  nRCell = 12;
  nRMilieu = 2;
  nRCellMilieu = 2;

  // Conversion coefficient between number of molecules and concentration in nM.
  nMConversion = 1.6605392767355127;
  nMConversionPower3 = pow(nMConversion,3);
}

//------------------------------------------------------------------------------

ChemicalLangevinComputeIncrementProfiling::~ChemicalLangevinComputeIncrementProfiling()
{}

//------------------------------------------------------------------------------

void ChemicalLangevinComputeIncrementProfiling::setFixedNCells(int nCells)
{
    cout << "NCells in the ChemicalLangevinComputeIncrementProfiling class = " << nCells << endl;
  nSCellstot = nCells*nSCell;
  nRtot = 2.*(nCells*(nRCell+nRCellMilieu) + nRMilieu);
  nStot = nCells*nSCell + nSMilieu;

  bNoise.resize(nStot,nRtot);
  bNoise = 0.;

  noiseIncrementVector.resize(nStot);
}


double ChemicalLangevinComputeIncrementProfiling::computeMean(Array<double,1>& array)
{
  return mean( array );
}

void ChemicalLangevinComputeIncrementProfiling::computeDeterministicIncrement
(
  Array<double,1>& xConcCells,
  Array<double,1>& xConcMilieu,
  const int nCells,
  const double volumeCell,
  const double volumeExt,
  const double chemicalLangevinTimeStep,
  double avgU,
  double avgV
)
{
  // Computing deterministic increments

  // Compute increments for cell species
  for (int iCell=0; iCell<nCells; iCell++)
  {
    xConcCells(iCell*nSCell + 0) += (0.5 + avgU*Input::globalParameter(0) - 1.*xConcCells(2*iCell) - 1.*Input::globalParameter(0)*xConcCells(2*iCell) + 1.25/(0.125 + xConcCells(1 + 2*iCell)*xConcCells(1 + 2*iCell)*xConcCells(1 + 2*iCell))) * chemicalLangevinTimeStep;
    xConcCells(iCell*nSCell + 1) += (0.5 + avgV*Input::globalParameter(0) + 1.25/(0.125 + xConcCells(2*iCell)*xConcCells(2*iCell)*xConcCells(2*iCell)) - 1.*xConcCells(1 + 2*iCell) - 1.*Input::globalParameter(0)*xConcCells(1 + 2*iCell)) * chemicalLangevinTimeStep;
  }

  // Compute increments for milieu species. We need to compute separately the increments
  // for the milieu reactions and for the cell-milieu reactions, which will be corrected for volume effects.

  // Compute increments for milieu species for milieu reactions.
  xConcMilieu(0) += (0.) * chemicalLangevinTimeStep;

  // Compute increments for milieu species for cell-milieu reactions.
  for (int iCell=0; iCell<nCells; iCell++)
  {
    xConcMilieu(0) += (0.) * chemicalLangevinTimeStep;
  }

}

void ChemicalLangevinComputeIncrementProfiling::generateRandomVector
(
  Array<double,1>& randomGaussianVector
)
{
  // Generating a vector of nRtot gaussian normal random numbers
  for (int i=0; i<nRtot; i++)
  {
    randomGaussianVector(i) = RandomNumberGenerator::getNormal();
  }
}


void ChemicalLangevinComputeIncrementProfiling::computeBNoise
(
  Array<double,2>& bNoise,
  Array<double,1>& xConcCells,
  Array<double,1>& xConcMilieu,
  double avgU,
  double avgV
)
{
  double noiseTermDAvgU = sqrt(avgU*Input::globalParameter(0));
  double noiseTermDAvgV = sqrt(avgV*Input::globalParameter(0));

//  bNoise.resize(nStot,nRtot);
//  bNoise = 0.;

  bNoise(0,0) = sqrt(0.5 + 1.25/(0.125 + xConcCells(1)*xConcCells(1)*xConcCells(1)));
  bNoise(0,2) = noiseTermDAvgU;
  bNoise(0,3) = -1.*sqrt(Input::globalParameter(0)*xConcCells(0));
  bNoise(0,701) = -1.*sqrt(xConcCells(0));
  bNoise(1,1) = sqrt(0.5 + 1.25/(0.125 + xConcCells(0)*xConcCells(0)*xConcCells(0)));
  bNoise(1,4) = noiseTermDAvgV;
  bNoise(1,5) = -1.*sqrt(Input::globalParameter(0)*xConcCells(1));
  bNoise(1,702) = -1.*sqrt(xConcCells(1));
  bNoise(2,6) = sqrt(0.5 + 1.25/(0.125 + xConcCells(3)*xConcCells(3)*xConcCells(3)));
  bNoise(2,8) = noiseTermDAvgU;
  bNoise(2,9) = -1.*sqrt(Input::globalParameter(0)*xConcCells(2));
  bNoise(2,707) = -1.*sqrt(xConcCells(2));
  bNoise(3,7) = sqrt(0.5 + 1.25/(0.125 + xConcCells(2)*xConcCells(2)*xConcCells(2)));
  bNoise(3,10) = noiseTermDAvgV;
  bNoise(3,11) = -1.*sqrt(Input::globalParameter(0)*xConcCells(3));
  bNoise(3,708) = -1.*sqrt(xConcCells(3));
  bNoise(4,12) = sqrt(0.5 + 1.25/(0.125 + xConcCells(5)*xConcCells(5)*xConcCells(5)));
  bNoise(4,14) = noiseTermDAvgU;
  bNoise(4,15) = -1.*sqrt(Input::globalParameter(0)*xConcCells(4));
  bNoise(4,713) = -1.*sqrt(xConcCells(4));
  bNoise(5,13) = sqrt(0.5 + 1.25/(0.125 + xConcCells(4)*xConcCells(4)*xConcCells(4)));
  bNoise(5,16) = noiseTermDAvgV;
  bNoise(5,17) = -1.*sqrt(Input::globalParameter(0)*xConcCells(5));
  bNoise(5,714) = -1.*sqrt(xConcCells(5));
  bNoise(6,18) = sqrt(0.5 + 1.25/(0.125 + xConcCells(7)*xConcCells(7)*xConcCells(7)));
  bNoise(6,20) = noiseTermDAvgU;
  bNoise(6,21) = -1.*sqrt(Input::globalParameter(0)*xConcCells(6));
  bNoise(6,719) = -1.*sqrt(xConcCells(6));
  bNoise(7,19) = sqrt(0.5 + 1.25/(0.125 + xConcCells(6)*xConcCells(6)*xConcCells(6)));
  bNoise(7,22) = noiseTermDAvgV;
  bNoise(7,23) = -1.*sqrt(Input::globalParameter(0)*xConcCells(7));
  bNoise(7,720) = -1.*sqrt(xConcCells(7));
  bNoise(8,24) = sqrt(0.5 + 1.25/(0.125 + xConcCells(9)*xConcCells(9)*xConcCells(9)));
  bNoise(8,26) = noiseTermDAvgU;
  bNoise(8,27) = -1.*sqrt(Input::globalParameter(0)*xConcCells(8));
  bNoise(8,725) = -1.*sqrt(xConcCells(8));
  bNoise(9,25) = sqrt(0.5 + 1.25/(0.125 + xConcCells(8)*xConcCells(8)*xConcCells(8)));
  bNoise(9,28) = noiseTermDAvgV;
  bNoise(9,29) = -1.*sqrt(Input::globalParameter(0)*xConcCells(9));
  bNoise(9,726) = -1.*sqrt(xConcCells(9));
  bNoise(10,30) = sqrt(0.5 + 1.25/(0.125 + xConcCells(11)*xConcCells(11)*xConcCells(11)));
  bNoise(10,32) = noiseTermDAvgU;
  bNoise(10,33) = -1.*sqrt(Input::globalParameter(0)*xConcCells(10));
  bNoise(10,731) = -1.*sqrt(xConcCells(10));
  bNoise(11,31) = sqrt(0.5 + 1.25/(0.125 + xConcCells(10)*xConcCells(10)*xConcCells(10)));
  bNoise(11,34) = noiseTermDAvgV;
  bNoise(11,35) = -1.*sqrt(Input::globalParameter(0)*xConcCells(11));
  bNoise(11,732) = -1.*sqrt(xConcCells(11));
  bNoise(12,36) = sqrt(0.5 + 1.25/(0.125 + xConcCells(13)*xConcCells(13)*xConcCells(13)));
  bNoise(12,38) = noiseTermDAvgU;
  bNoise(12,39) = -1.*sqrt(Input::globalParameter(0)*xConcCells(12));
  bNoise(12,737) = -1.*sqrt(xConcCells(12));
  bNoise(13,37) = sqrt(0.5 + 1.25/(0.125 + xConcCells(12)*xConcCells(12)*xConcCells(12)));
  bNoise(13,40) = noiseTermDAvgV;
  bNoise(13,41) = -1.*sqrt(Input::globalParameter(0)*xConcCells(13));
  bNoise(13,738) = -1.*sqrt(xConcCells(13));
  bNoise(14,42) = sqrt(0.5 + 1.25/(0.125 + xConcCells(15)*xConcCells(15)*xConcCells(15)));
  bNoise(14,44) = noiseTermDAvgU;
  bNoise(14,45) = -1.*sqrt(Input::globalParameter(0)*xConcCells(14));
  bNoise(14,743) = -1.*sqrt(xConcCells(14));
  bNoise(15,43) = sqrt(0.5 + 1.25/(0.125 + xConcCells(14)*xConcCells(14)*xConcCells(14)));
  bNoise(15,46) = noiseTermDAvgV;
  bNoise(15,47) = -1.*sqrt(Input::globalParameter(0)*xConcCells(15));
  bNoise(15,744) = -1.*sqrt(xConcCells(15));
  bNoise(16,48) = sqrt(0.5 + 1.25/(0.125 + xConcCells(17)*xConcCells(17)*xConcCells(17)));
  bNoise(16,50) = noiseTermDAvgU;
  bNoise(16,51) = -1.*sqrt(Input::globalParameter(0)*xConcCells(16));
  bNoise(16,749) = -1.*sqrt(xConcCells(16));
  bNoise(17,49) = sqrt(0.5 + 1.25/(0.125 + xConcCells(16)*xConcCells(16)*xConcCells(16)));
  bNoise(17,52) = noiseTermDAvgV;
  bNoise(17,53) = -1.*sqrt(Input::globalParameter(0)*xConcCells(17));
  bNoise(17,750) = -1.*sqrt(xConcCells(17));
  bNoise(18,54) = sqrt(0.5 + 1.25/(0.125 + xConcCells(19)*xConcCells(19)*xConcCells(19)));
  bNoise(18,56) = noiseTermDAvgU;
  bNoise(18,57) = -1.*sqrt(Input::globalParameter(0)*xConcCells(18));
  bNoise(18,755) = -1.*sqrt(xConcCells(18));
  bNoise(19,55) = sqrt(0.5 + 1.25/(0.125 + xConcCells(18)*xConcCells(18)*xConcCells(18)));
  bNoise(19,58) = noiseTermDAvgV;
  bNoise(19,59) = -1.*sqrt(Input::globalParameter(0)*xConcCells(19));
  bNoise(19,756) = -1.*sqrt(xConcCells(19));
  bNoise(20,60) = sqrt(0.5 + 1.25/(0.125 + xConcCells(21)*xConcCells(21)*xConcCells(21)));
  bNoise(20,62) = noiseTermDAvgU;
  bNoise(20,63) = -1.*sqrt(Input::globalParameter(0)*xConcCells(20));
  bNoise(20,761) = -1.*sqrt(xConcCells(20));
  bNoise(21,61) = sqrt(0.5 + 1.25/(0.125 + xConcCells(20)*xConcCells(20)*xConcCells(20)));
  bNoise(21,64) = noiseTermDAvgV;
  bNoise(21,65) = -1.*sqrt(Input::globalParameter(0)*xConcCells(21));
  bNoise(21,762) = -1.*sqrt(xConcCells(21));
  bNoise(22,66) = sqrt(0.5 + 1.25/(0.125 + xConcCells(23)*xConcCells(23)*xConcCells(23)));
  bNoise(22,68) = noiseTermDAvgU;
  bNoise(22,69) = -1.*sqrt(Input::globalParameter(0)*xConcCells(22));
  bNoise(22,767) = -1.*sqrt(xConcCells(22));
  bNoise(23,67) = sqrt(0.5 + 1.25/(0.125 + xConcCells(22)*xConcCells(22)*xConcCells(22)));
  bNoise(23,70) = noiseTermDAvgV;
  bNoise(23,71) = -1.*sqrt(Input::globalParameter(0)*xConcCells(23));
  bNoise(23,768) = -1.*sqrt(xConcCells(23));
  bNoise(24,72) = sqrt(0.5 + 1.25/(0.125 + xConcCells(25)*xConcCells(25)*xConcCells(25)));
  bNoise(24,74) = noiseTermDAvgU;
  bNoise(24,75) = -1.*sqrt(Input::globalParameter(0)*xConcCells(24));
  bNoise(24,773) = -1.*sqrt(xConcCells(24));
  bNoise(25,73) = sqrt(0.5 + 1.25/(0.125 + xConcCells(24)*xConcCells(24)*xConcCells(24)));
  bNoise(25,76) = noiseTermDAvgV;
  bNoise(25,77) = -1.*sqrt(Input::globalParameter(0)*xConcCells(25));
  bNoise(25,774) = -1.*sqrt(xConcCells(25));
  bNoise(26,78) = sqrt(0.5 + 1.25/(0.125 + xConcCells(27)*xConcCells(27)*xConcCells(27)));
  bNoise(26,80) = noiseTermDAvgU;
  bNoise(26,81) = -1.*sqrt(Input::globalParameter(0)*xConcCells(26));
  bNoise(26,779) = -1.*sqrt(xConcCells(26));
  bNoise(27,79) = sqrt(0.5 + 1.25/(0.125 + xConcCells(26)*xConcCells(26)*xConcCells(26)));
  bNoise(27,82) = noiseTermDAvgV;
  bNoise(27,83) = -1.*sqrt(Input::globalParameter(0)*xConcCells(27));
  bNoise(27,780) = -1.*sqrt(xConcCells(27));
  bNoise(28,84) = sqrt(0.5 + 1.25/(0.125 + xConcCells(29)*xConcCells(29)*xConcCells(29)));
  bNoise(28,86) = noiseTermDAvgU;
  bNoise(28,87) = -1.*sqrt(Input::globalParameter(0)*xConcCells(28));
  bNoise(28,785) = -1.*sqrt(xConcCells(28));
  bNoise(29,85) = sqrt(0.5 + 1.25/(0.125 + xConcCells(28)*xConcCells(28)*xConcCells(28)));
  bNoise(29,88) = noiseTermDAvgV;
  bNoise(29,89) = -1.*sqrt(Input::globalParameter(0)*xConcCells(29));
  bNoise(29,786) = -1.*sqrt(xConcCells(29));
  bNoise(30,90) = sqrt(0.5 + 1.25/(0.125 + xConcCells(31)*xConcCells(31)*xConcCells(31)));
  bNoise(30,92) = noiseTermDAvgU;
  bNoise(30,93) = -1.*sqrt(Input::globalParameter(0)*xConcCells(30));
  bNoise(30,791) = -1.*sqrt(xConcCells(30));
  bNoise(31,91) = sqrt(0.5 + 1.25/(0.125 + xConcCells(30)*xConcCells(30)*xConcCells(30)));
  bNoise(31,94) = noiseTermDAvgV;
  bNoise(31,95) = -1.*sqrt(Input::globalParameter(0)*xConcCells(31));
  bNoise(31,792) = -1.*sqrt(xConcCells(31));
  bNoise(32,96) = sqrt(0.5 + 1.25/(0.125 + xConcCells(33)*xConcCells(33)*xConcCells(33)));
  bNoise(32,98) = noiseTermDAvgU;
  bNoise(32,99) = -1.*sqrt(Input::globalParameter(0)*xConcCells(32));
  bNoise(32,797) = -1.*sqrt(xConcCells(32));
  bNoise(33,97) = sqrt(0.5 + 1.25/(0.125 + xConcCells(32)*xConcCells(32)*xConcCells(32)));
  bNoise(33,100) = noiseTermDAvgV;
  bNoise(33,101) = -1.*sqrt(Input::globalParameter(0)*xConcCells(33));
  bNoise(33,798) = -1.*sqrt(xConcCells(33));
  bNoise(34,102) = sqrt(0.5 + 1.25/(0.125 + xConcCells(35)*xConcCells(35)*xConcCells(35)));
  bNoise(34,104) = noiseTermDAvgU;
  bNoise(34,105) = -1.*sqrt(Input::globalParameter(0)*xConcCells(34));
  bNoise(34,803) = -1.*sqrt(xConcCells(34));
  bNoise(35,103) = sqrt(0.5 + 1.25/(0.125 + xConcCells(34)*xConcCells(34)*xConcCells(34)));
  bNoise(35,106) = noiseTermDAvgV;
  bNoise(35,107) = -1.*sqrt(Input::globalParameter(0)*xConcCells(35));
  bNoise(35,804) = -1.*sqrt(xConcCells(35));
  bNoise(36,108) = sqrt(0.5 + 1.25/(0.125 + xConcCells(37)*xConcCells(37)*xConcCells(37)));
  bNoise(36,110) = noiseTermDAvgU;
  bNoise(36,111) = -1.*sqrt(Input::globalParameter(0)*xConcCells(36));
  bNoise(36,809) = -1.*sqrt(xConcCells(36));
  bNoise(37,109) = sqrt(0.5 + 1.25/(0.125 + xConcCells(36)*xConcCells(36)*xConcCells(36)));
  bNoise(37,112) = noiseTermDAvgV;
  bNoise(37,113) = -1.*sqrt(Input::globalParameter(0)*xConcCells(37));
  bNoise(37,810) = -1.*sqrt(xConcCells(37));
  bNoise(38,114) = sqrt(0.5 + 1.25/(0.125 + xConcCells(39)*xConcCells(39)*xConcCells(39)));
  bNoise(38,116) = noiseTermDAvgU;
  bNoise(38,117) = -1.*sqrt(Input::globalParameter(0)*xConcCells(38));
  bNoise(38,815) = -1.*sqrt(xConcCells(38));
  bNoise(39,115) = sqrt(0.5 + 1.25/(0.125 + xConcCells(38)*xConcCells(38)*xConcCells(38)));
  bNoise(39,118) = noiseTermDAvgV;
  bNoise(39,119) = -1.*sqrt(Input::globalParameter(0)*xConcCells(39));
  bNoise(39,816) = -1.*sqrt(xConcCells(39));
  bNoise(40,120) = sqrt(0.5 + 1.25/(0.125 + xConcCells(41)*xConcCells(41)*xConcCells(41)));
  bNoise(40,122) = noiseTermDAvgU;
  bNoise(40,123) = -1.*sqrt(Input::globalParameter(0)*xConcCells(40));
  bNoise(40,821) = -1.*sqrt(xConcCells(40));
  bNoise(41,121) = sqrt(0.5 + 1.25/(0.125 + xConcCells(40)*xConcCells(40)*xConcCells(40)));
  bNoise(41,124) = noiseTermDAvgV;
  bNoise(41,125) = -1.*sqrt(Input::globalParameter(0)*xConcCells(41));
  bNoise(41,822) = -1.*sqrt(xConcCells(41));
  bNoise(42,126) = sqrt(0.5 + 1.25/(0.125 + xConcCells(43)*xConcCells(43)*xConcCells(43)));
  bNoise(42,128) = noiseTermDAvgU;
  bNoise(42,129) = -1.*sqrt(Input::globalParameter(0)*xConcCells(42));
  bNoise(42,827) = -1.*sqrt(xConcCells(42));
  bNoise(43,127) = sqrt(0.5 + 1.25/(0.125 + xConcCells(42)*xConcCells(42)*xConcCells(42)));
  bNoise(43,130) = noiseTermDAvgV;
  bNoise(43,131) = -1.*sqrt(Input::globalParameter(0)*xConcCells(43));
  bNoise(43,828) = -1.*sqrt(xConcCells(43));
  bNoise(44,132) = sqrt(0.5 + 1.25/(0.125 + xConcCells(45)*xConcCells(45)*xConcCells(45)));
  bNoise(44,134) = noiseTermDAvgU;
  bNoise(44,135) = -1.*sqrt(Input::globalParameter(0)*xConcCells(44));
  bNoise(44,833) = -1.*sqrt(xConcCells(44));
  bNoise(45,133) = sqrt(0.5 + 1.25/(0.125 + xConcCells(44)*xConcCells(44)*xConcCells(44)));
  bNoise(45,136) = noiseTermDAvgV;
  bNoise(45,137) = -1.*sqrt(Input::globalParameter(0)*xConcCells(45));
  bNoise(45,834) = -1.*sqrt(xConcCells(45));
  bNoise(46,138) = sqrt(0.5 + 1.25/(0.125 + xConcCells(47)*xConcCells(47)*xConcCells(47)));
  bNoise(46,140) = noiseTermDAvgU;
  bNoise(46,141) = -1.*sqrt(Input::globalParameter(0)*xConcCells(46));
  bNoise(46,839) = -1.*sqrt(xConcCells(46));
  bNoise(47,139) = sqrt(0.5 + 1.25/(0.125 + xConcCells(46)*xConcCells(46)*xConcCells(46)));
  bNoise(47,142) = noiseTermDAvgV;
  bNoise(47,143) = -1.*sqrt(Input::globalParameter(0)*xConcCells(47));
  bNoise(47,840) = -1.*sqrt(xConcCells(47));
  bNoise(48,144) = sqrt(0.5 + 1.25/(0.125 + xConcCells(49)*xConcCells(49)*xConcCells(49)));
  bNoise(48,146) = noiseTermDAvgU;
  bNoise(48,147) = -1.*sqrt(Input::globalParameter(0)*xConcCells(48));
  bNoise(48,845) = -1.*sqrt(xConcCells(48));
  bNoise(49,145) = sqrt(0.5 + 1.25/(0.125 + xConcCells(48)*xConcCells(48)*xConcCells(48)));
  bNoise(49,148) = noiseTermDAvgV;
  bNoise(49,149) = -1.*sqrt(Input::globalParameter(0)*xConcCells(49));
  bNoise(49,846) = -1.*sqrt(xConcCells(49));
  bNoise(50,150) = sqrt(0.5 + 1.25/(0.125 + xConcCells(51)*xConcCells(51)*xConcCells(51)));
  bNoise(50,152) = noiseTermDAvgU;
  bNoise(50,153) = -1.*sqrt(Input::globalParameter(0)*xConcCells(50));
  bNoise(50,851) = -1.*sqrt(xConcCells(50));
  bNoise(51,151) = sqrt(0.5 + 1.25/(0.125 + xConcCells(50)*xConcCells(50)*xConcCells(50)));
  bNoise(51,154) = noiseTermDAvgV;
  bNoise(51,155) = -1.*sqrt(Input::globalParameter(0)*xConcCells(51));
  bNoise(51,852) = -1.*sqrt(xConcCells(51));
  bNoise(52,156) = sqrt(0.5 + 1.25/(0.125 + xConcCells(53)*xConcCells(53)*xConcCells(53)));
  bNoise(52,158) = noiseTermDAvgU;
  bNoise(52,159) = -1.*sqrt(Input::globalParameter(0)*xConcCells(52));
  bNoise(52,857) = -1.*sqrt(xConcCells(52));
  bNoise(53,157) = sqrt(0.5 + 1.25/(0.125 + xConcCells(52)*xConcCells(52)*xConcCells(52)));
  bNoise(53,160) = noiseTermDAvgV;
  bNoise(53,161) = -1.*sqrt(Input::globalParameter(0)*xConcCells(53));
  bNoise(53,858) = -1.*sqrt(xConcCells(53));
  bNoise(54,162) = sqrt(0.5 + 1.25/(0.125 + xConcCells(55)*xConcCells(55)*xConcCells(55)));
  bNoise(54,164) = noiseTermDAvgU;
  bNoise(54,165) = -1.*sqrt(Input::globalParameter(0)*xConcCells(54));
  bNoise(54,863) = -1.*sqrt(xConcCells(54));
  bNoise(55,163) = sqrt(0.5 + 1.25/(0.125 + xConcCells(54)*xConcCells(54)*xConcCells(54)));
  bNoise(55,166) = noiseTermDAvgV;
  bNoise(55,167) = -1.*sqrt(Input::globalParameter(0)*xConcCells(55));
  bNoise(55,864) = -1.*sqrt(xConcCells(55));
  bNoise(56,168) = sqrt(0.5 + 1.25/(0.125 + xConcCells(57)*xConcCells(57)*xConcCells(57)));
  bNoise(56,170) = noiseTermDAvgU;
  bNoise(56,171) = -1.*sqrt(Input::globalParameter(0)*xConcCells(56));
  bNoise(56,869) = -1.*sqrt(xConcCells(56));
  bNoise(57,169) = sqrt(0.5 + 1.25/(0.125 + xConcCells(56)*xConcCells(56)*xConcCells(56)));
  bNoise(57,172) = noiseTermDAvgV;
  bNoise(57,173) = -1.*sqrt(Input::globalParameter(0)*xConcCells(57));
  bNoise(57,870) = -1.*sqrt(xConcCells(57));
  bNoise(58,174) = sqrt(0.5 + 1.25/(0.125 + xConcCells(59)*xConcCells(59)*xConcCells(59)));
  bNoise(58,176) = noiseTermDAvgU;
  bNoise(58,177) = -1.*sqrt(Input::globalParameter(0)*xConcCells(58));
  bNoise(58,875) = -1.*sqrt(xConcCells(58));
  bNoise(59,175) = sqrt(0.5 + 1.25/(0.125 + xConcCells(58)*xConcCells(58)*xConcCells(58)));
  bNoise(59,178) = noiseTermDAvgV;
  bNoise(59,179) = -1.*sqrt(Input::globalParameter(0)*xConcCells(59));
  bNoise(59,876) = -1.*sqrt(xConcCells(59));
  bNoise(60,180) = sqrt(0.5 + 1.25/(0.125 + xConcCells(61)*xConcCells(61)*xConcCells(61)));
  bNoise(60,182) = noiseTermDAvgU;
  bNoise(60,183) = -1.*sqrt(Input::globalParameter(0)*xConcCells(60));
  bNoise(60,881) = -1.*sqrt(xConcCells(60));
  bNoise(61,181) = sqrt(0.5 + 1.25/(0.125 + xConcCells(60)*xConcCells(60)*xConcCells(60)));
  bNoise(61,184) = noiseTermDAvgV;
  bNoise(61,185) = -1.*sqrt(Input::globalParameter(0)*xConcCells(61));
  bNoise(61,882) = -1.*sqrt(xConcCells(61));
  bNoise(62,186) = sqrt(0.5 + 1.25/(0.125 + xConcCells(63)*xConcCells(63)*xConcCells(63)));
  bNoise(62,188) = noiseTermDAvgU;
  bNoise(62,189) = -1.*sqrt(Input::globalParameter(0)*xConcCells(62));
  bNoise(62,887) = -1.*sqrt(xConcCells(62));
  bNoise(63,187) = sqrt(0.5 + 1.25/(0.125 + xConcCells(62)*xConcCells(62)*xConcCells(62)));
  bNoise(63,190) = noiseTermDAvgV;
  bNoise(63,191) = -1.*sqrt(Input::globalParameter(0)*xConcCells(63));
  bNoise(63,888) = -1.*sqrt(xConcCells(63));
  bNoise(64,192) = sqrt(0.5 + 1.25/(0.125 + xConcCells(65)*xConcCells(65)*xConcCells(65)));
  bNoise(64,194) = noiseTermDAvgU;
  bNoise(64,195) = -1.*sqrt(Input::globalParameter(0)*xConcCells(64));
  bNoise(64,893) = -1.*sqrt(xConcCells(64));
  bNoise(65,193) = sqrt(0.5 + 1.25/(0.125 + xConcCells(64)*xConcCells(64)*xConcCells(64)));
  bNoise(65,196) = noiseTermDAvgV;
  bNoise(65,197) = -1.*sqrt(Input::globalParameter(0)*xConcCells(65));
  bNoise(65,894) = -1.*sqrt(xConcCells(65));
  bNoise(66,198) = sqrt(0.5 + 1.25/(0.125 + xConcCells(67)*xConcCells(67)*xConcCells(67)));
  bNoise(66,200) = noiseTermDAvgU;
  bNoise(66,201) = -1.*sqrt(Input::globalParameter(0)*xConcCells(66));
  bNoise(66,899) = -1.*sqrt(xConcCells(66));
  bNoise(67,199) = sqrt(0.5 + 1.25/(0.125 + xConcCells(66)*xConcCells(66)*xConcCells(66)));
  bNoise(67,202) = noiseTermDAvgV;
  bNoise(67,203) = -1.*sqrt(Input::globalParameter(0)*xConcCells(67));
  bNoise(67,900) = -1.*sqrt(xConcCells(67));
  bNoise(68,204) = sqrt(0.5 + 1.25/(0.125 + xConcCells(69)*xConcCells(69)*xConcCells(69)));
  bNoise(68,206) = noiseTermDAvgU;
  bNoise(68,207) = -1.*sqrt(Input::globalParameter(0)*xConcCells(68));
  bNoise(68,905) = -1.*sqrt(xConcCells(68));
  bNoise(69,205) = sqrt(0.5 + 1.25/(0.125 + xConcCells(68)*xConcCells(68)*xConcCells(68)));
  bNoise(69,208) = noiseTermDAvgV;
  bNoise(69,209) = -1.*sqrt(Input::globalParameter(0)*xConcCells(69));
  bNoise(69,906) = -1.*sqrt(xConcCells(69));
  bNoise(70,210) = sqrt(0.5 + 1.25/(0.125 + xConcCells(71)*xConcCells(71)*xConcCells(71)));
  bNoise(70,212) = noiseTermDAvgU;
  bNoise(70,213) = -1.*sqrt(Input::globalParameter(0)*xConcCells(70));
  bNoise(70,911) = -1.*sqrt(xConcCells(70));
  bNoise(71,211) = sqrt(0.5 + 1.25/(0.125 + xConcCells(70)*xConcCells(70)*xConcCells(70)));
  bNoise(71,214) = noiseTermDAvgV;
  bNoise(71,215) = -1.*sqrt(Input::globalParameter(0)*xConcCells(71));
  bNoise(71,912) = -1.*sqrt(xConcCells(71));
  bNoise(72,216) = sqrt(0.5 + 1.25/(0.125 + xConcCells(73)*xConcCells(73)*xConcCells(73)));
  bNoise(72,218) = noiseTermDAvgU;
  bNoise(72,219) = -1.*sqrt(Input::globalParameter(0)*xConcCells(72));
  bNoise(72,917) = -1.*sqrt(xConcCells(72));
  bNoise(73,217) = sqrt(0.5 + 1.25/(0.125 + xConcCells(72)*xConcCells(72)*xConcCells(72)));
  bNoise(73,220) = noiseTermDAvgV;
  bNoise(73,221) = -1.*sqrt(Input::globalParameter(0)*xConcCells(73));
  bNoise(73,918) = -1.*sqrt(xConcCells(73));
  bNoise(74,222) = sqrt(0.5 + 1.25/(0.125 + xConcCells(75)*xConcCells(75)*xConcCells(75)));
  bNoise(74,224) = noiseTermDAvgU;
  bNoise(74,225) = -1.*sqrt(Input::globalParameter(0)*xConcCells(74));
  bNoise(74,923) = -1.*sqrt(xConcCells(74));
  bNoise(75,223) = sqrt(0.5 + 1.25/(0.125 + xConcCells(74)*xConcCells(74)*xConcCells(74)));
  bNoise(75,226) = noiseTermDAvgV;
  bNoise(75,227) = -1.*sqrt(Input::globalParameter(0)*xConcCells(75));
  bNoise(75,924) = -1.*sqrt(xConcCells(75));
  bNoise(76,228) = sqrt(0.5 + 1.25/(0.125 + xConcCells(77)*xConcCells(77)*xConcCells(77)));
  bNoise(76,230) = noiseTermDAvgU;
  bNoise(76,231) = -1.*sqrt(Input::globalParameter(0)*xConcCells(76));
  bNoise(76,929) = -1.*sqrt(xConcCells(76));
  bNoise(77,229) = sqrt(0.5 + 1.25/(0.125 + xConcCells(76)*xConcCells(76)*xConcCells(76)));
  bNoise(77,232) = noiseTermDAvgV;
  bNoise(77,233) = -1.*sqrt(Input::globalParameter(0)*xConcCells(77));
  bNoise(77,930) = -1.*sqrt(xConcCells(77));
  bNoise(78,234) = sqrt(0.5 + 1.25/(0.125 + xConcCells(79)*xConcCells(79)*xConcCells(79)));
  bNoise(78,236) = noiseTermDAvgU;
  bNoise(78,237) = -1.*sqrt(Input::globalParameter(0)*xConcCells(78));
  bNoise(78,935) = -1.*sqrt(xConcCells(78));
  bNoise(79,235) = sqrt(0.5 + 1.25/(0.125 + xConcCells(78)*xConcCells(78)*xConcCells(78)));
  bNoise(79,238) = noiseTermDAvgV;
  bNoise(79,239) = -1.*sqrt(Input::globalParameter(0)*xConcCells(79));
  bNoise(79,936) = -1.*sqrt(xConcCells(79));
  bNoise(80,240) = sqrt(0.5 + 1.25/(0.125 + xConcCells(81)*xConcCells(81)*xConcCells(81)));
  bNoise(80,242) = noiseTermDAvgU;
  bNoise(80,243) = -1.*sqrt(Input::globalParameter(0)*xConcCells(80));
  bNoise(80,941) = -1.*sqrt(xConcCells(80));
  bNoise(81,241) = sqrt(0.5 + 1.25/(0.125 + xConcCells(80)*xConcCells(80)*xConcCells(80)));
  bNoise(81,244) = noiseTermDAvgV;
  bNoise(81,245) = -1.*sqrt(Input::globalParameter(0)*xConcCells(81));
  bNoise(81,942) = -1.*sqrt(xConcCells(81));
  bNoise(82,246) = sqrt(0.5 + 1.25/(0.125 + xConcCells(83)*xConcCells(83)*xConcCells(83)));
  bNoise(82,248) = noiseTermDAvgU;
  bNoise(82,249) = -1.*sqrt(Input::globalParameter(0)*xConcCells(82));
  bNoise(82,947) = -1.*sqrt(xConcCells(82));
  bNoise(83,247) = sqrt(0.5 + 1.25/(0.125 + xConcCells(82)*xConcCells(82)*xConcCells(82)));
  bNoise(83,250) = noiseTermDAvgV;
  bNoise(83,251) = -1.*sqrt(Input::globalParameter(0)*xConcCells(83));
  bNoise(83,948) = -1.*sqrt(xConcCells(83));
  bNoise(84,252) = sqrt(0.5 + 1.25/(0.125 + xConcCells(85)*xConcCells(85)*xConcCells(85)));
  bNoise(84,254) = noiseTermDAvgU;
  bNoise(84,255) = -1.*sqrt(Input::globalParameter(0)*xConcCells(84));
  bNoise(84,953) = -1.*sqrt(xConcCells(84));
  bNoise(85,253) = sqrt(0.5 + 1.25/(0.125 + xConcCells(84)*xConcCells(84)*xConcCells(84)));
  bNoise(85,256) = noiseTermDAvgV;
  bNoise(85,257) = -1.*sqrt(Input::globalParameter(0)*xConcCells(85));
  bNoise(85,954) = -1.*sqrt(xConcCells(85));
  bNoise(86,258) = sqrt(0.5 + 1.25/(0.125 + xConcCells(87)*xConcCells(87)*xConcCells(87)));
  bNoise(86,260) = noiseTermDAvgU;
  bNoise(86,261) = -1.*sqrt(Input::globalParameter(0)*xConcCells(86));
  bNoise(86,959) = -1.*sqrt(xConcCells(86));
  bNoise(87,259) = sqrt(0.5 + 1.25/(0.125 + xConcCells(86)*xConcCells(86)*xConcCells(86)));
  bNoise(87,262) = noiseTermDAvgV;
  bNoise(87,263) = -1.*sqrt(Input::globalParameter(0)*xConcCells(87));
  bNoise(87,960) = -1.*sqrt(xConcCells(87));
  bNoise(88,264) = sqrt(0.5 + 1.25/(0.125 + xConcCells(89)*xConcCells(89)*xConcCells(89)));
  bNoise(88,266) = noiseTermDAvgU;
  bNoise(88,267) = -1.*sqrt(Input::globalParameter(0)*xConcCells(88));
  bNoise(88,965) = -1.*sqrt(xConcCells(88));
  bNoise(89,265) = sqrt(0.5 + 1.25/(0.125 + xConcCells(88)*xConcCells(88)*xConcCells(88)));
  bNoise(89,268) = noiseTermDAvgV;
  bNoise(89,269) = -1.*sqrt(Input::globalParameter(0)*xConcCells(89));
  bNoise(89,966) = -1.*sqrt(xConcCells(89));
  bNoise(90,270) = sqrt(0.5 + 1.25/(0.125 + xConcCells(91)*xConcCells(91)*xConcCells(91)));
  bNoise(90,272) = noiseTermDAvgU;
  bNoise(90,273) = -1.*sqrt(Input::globalParameter(0)*xConcCells(90));
  bNoise(90,971) = -1.*sqrt(xConcCells(90));
  bNoise(91,271) = sqrt(0.5 + 1.25/(0.125 + xConcCells(90)*xConcCells(90)*xConcCells(90)));
  bNoise(91,274) = noiseTermDAvgV;
  bNoise(91,275) = -1.*sqrt(Input::globalParameter(0)*xConcCells(91));
  bNoise(91,972) = -1.*sqrt(xConcCells(91));
  bNoise(92,276) = sqrt(0.5 + 1.25/(0.125 + xConcCells(93)*xConcCells(93)*xConcCells(93)));
  bNoise(92,278) = noiseTermDAvgU;
  bNoise(92,279) = -1.*sqrt(Input::globalParameter(0)*xConcCells(92));
  bNoise(92,977) = -1.*sqrt(xConcCells(92));
  bNoise(93,277) = sqrt(0.5 + 1.25/(0.125 + xConcCells(92)*xConcCells(92)*xConcCells(92)));
  bNoise(93,280) = noiseTermDAvgV;
  bNoise(93,281) = -1.*sqrt(Input::globalParameter(0)*xConcCells(93));
  bNoise(93,978) = -1.*sqrt(xConcCells(93));
  bNoise(94,282) = sqrt(0.5 + 1.25/(0.125 + xConcCells(95)*xConcCells(95)*xConcCells(95)));
  bNoise(94,284) = noiseTermDAvgU;
  bNoise(94,285) = -1.*sqrt(Input::globalParameter(0)*xConcCells(94));
  bNoise(94,983) = -1.*sqrt(xConcCells(94));
  bNoise(95,283) = sqrt(0.5 + 1.25/(0.125 + xConcCells(94)*xConcCells(94)*xConcCells(94)));
  bNoise(95,286) = noiseTermDAvgV;
  bNoise(95,287) = -1.*sqrt(Input::globalParameter(0)*xConcCells(95));
  bNoise(95,984) = -1.*sqrt(xConcCells(95));
  bNoise(96,288) = sqrt(0.5 + 1.25/(0.125 + xConcCells(97)*xConcCells(97)*xConcCells(97)));
  bNoise(96,290) = noiseTermDAvgU;
  bNoise(96,291) = -1.*sqrt(Input::globalParameter(0)*xConcCells(96));
  bNoise(96,989) = -1.*sqrt(xConcCells(96));
  bNoise(97,289) = sqrt(0.5 + 1.25/(0.125 + xConcCells(96)*xConcCells(96)*xConcCells(96)));
  bNoise(97,292) = noiseTermDAvgV;
  bNoise(97,293) = -1.*sqrt(Input::globalParameter(0)*xConcCells(97));
  bNoise(97,990) = -1.*sqrt(xConcCells(97));
  bNoise(98,294) = sqrt(0.5 + 1.25/(0.125 + xConcCells(99)*xConcCells(99)*xConcCells(99)));
  bNoise(98,296) = noiseTermDAvgU;
  bNoise(98,297) = -1.*sqrt(Input::globalParameter(0)*xConcCells(98));
  bNoise(98,995) = -1.*sqrt(xConcCells(98));
  bNoise(99,295) = sqrt(0.5 + 1.25/(0.125 + xConcCells(98)*xConcCells(98)*xConcCells(98)));
  bNoise(99,298) = noiseTermDAvgV;
  bNoise(99,299) = -1.*sqrt(Input::globalParameter(0)*xConcCells(99));
  bNoise(99,996) = -1.*sqrt(xConcCells(99));
  bNoise(100,300) = sqrt(0.5 + 1.25/(0.125 + xConcCells(101)*xConcCells(101)*xConcCells(101)));
  bNoise(100,302) = noiseTermDAvgU;
  bNoise(100,303) = -1.*sqrt(Input::globalParameter(0)*xConcCells(100));
  bNoise(100,1001) = -1.*sqrt(xConcCells(100));
  bNoise(101,301) = sqrt(0.5 + 1.25/(0.125 + xConcCells(100)*xConcCells(100)*xConcCells(100)));
  bNoise(101,304) = noiseTermDAvgV;
  bNoise(101,305) = -1.*sqrt(Input::globalParameter(0)*xConcCells(101));
  bNoise(101,1002) = -1.*sqrt(xConcCells(101));
  bNoise(102,306) = sqrt(0.5 + 1.25/(0.125 + xConcCells(103)*xConcCells(103)*xConcCells(103)));
  bNoise(102,308) = noiseTermDAvgU;
  bNoise(102,309) = -1.*sqrt(Input::globalParameter(0)*xConcCells(102));
  bNoise(102,1007) = -1.*sqrt(xConcCells(102));
  bNoise(103,307) = sqrt(0.5 + 1.25/(0.125 + xConcCells(102)*xConcCells(102)*xConcCells(102)));
  bNoise(103,310) = noiseTermDAvgV;
  bNoise(103,311) = -1.*sqrt(Input::globalParameter(0)*xConcCells(103));
  bNoise(103,1008) = -1.*sqrt(xConcCells(103));
  bNoise(104,312) = sqrt(0.5 + 1.25/(0.125 + xConcCells(105)*xConcCells(105)*xConcCells(105)));
  bNoise(104,314) = noiseTermDAvgU;
  bNoise(104,315) = -1.*sqrt(Input::globalParameter(0)*xConcCells(104));
  bNoise(104,1013) = -1.*sqrt(xConcCells(104));
  bNoise(105,313) = sqrt(0.5 + 1.25/(0.125 + xConcCells(104)*xConcCells(104)*xConcCells(104)));
  bNoise(105,316) = noiseTermDAvgV;
  bNoise(105,317) = -1.*sqrt(Input::globalParameter(0)*xConcCells(105));
  bNoise(105,1014) = -1.*sqrt(xConcCells(105));
  bNoise(106,318) = sqrt(0.5 + 1.25/(0.125 + xConcCells(107)*xConcCells(107)*xConcCells(107)));
  bNoise(106,320) = noiseTermDAvgU;
  bNoise(106,321) = -1.*sqrt(Input::globalParameter(0)*xConcCells(106));
  bNoise(106,1019) = -1.*sqrt(xConcCells(106));
  bNoise(107,319) = sqrt(0.5 + 1.25/(0.125 + xConcCells(106)*xConcCells(106)*xConcCells(106)));
  bNoise(107,322) = noiseTermDAvgV;
  bNoise(107,323) = -1.*sqrt(Input::globalParameter(0)*xConcCells(107));
  bNoise(107,1020) = -1.*sqrt(xConcCells(107));
  bNoise(108,324) = sqrt(0.5 + 1.25/(0.125 + xConcCells(109)*xConcCells(109)*xConcCells(109)));
  bNoise(108,326) = noiseTermDAvgU;
  bNoise(108,327) = -1.*sqrt(Input::globalParameter(0)*xConcCells(108));
  bNoise(108,1025) = -1.*sqrt(xConcCells(108));
  bNoise(109,325) = sqrt(0.5 + 1.25/(0.125 + xConcCells(108)*xConcCells(108)*xConcCells(108)));
  bNoise(109,328) = noiseTermDAvgV;
  bNoise(109,329) = -1.*sqrt(Input::globalParameter(0)*xConcCells(109));
  bNoise(109,1026) = -1.*sqrt(xConcCells(109));
  bNoise(110,330) = sqrt(0.5 + 1.25/(0.125 + xConcCells(111)*xConcCells(111)*xConcCells(111)));
  bNoise(110,332) = noiseTermDAvgU;
  bNoise(110,333) = -1.*sqrt(Input::globalParameter(0)*xConcCells(110));
  bNoise(110,1031) = -1.*sqrt(xConcCells(110));
  bNoise(111,331) = sqrt(0.5 + 1.25/(0.125 + xConcCells(110)*xConcCells(110)*xConcCells(110)));
  bNoise(111,334) = noiseTermDAvgV;
  bNoise(111,335) = -1.*sqrt(Input::globalParameter(0)*xConcCells(111));
  bNoise(111,1032) = -1.*sqrt(xConcCells(111));
  bNoise(112,336) = sqrt(0.5 + 1.25/(0.125 + xConcCells(113)*xConcCells(113)*xConcCells(113)));
  bNoise(112,338) = noiseTermDAvgU;
  bNoise(112,339) = -1.*sqrt(Input::globalParameter(0)*xConcCells(112));
  bNoise(112,1037) = -1.*sqrt(xConcCells(112));
  bNoise(113,337) = sqrt(0.5 + 1.25/(0.125 + xConcCells(112)*xConcCells(112)*xConcCells(112)));
  bNoise(113,340) = noiseTermDAvgV;
  bNoise(113,341) = -1.*sqrt(Input::globalParameter(0)*xConcCells(113));
  bNoise(113,1038) = -1.*sqrt(xConcCells(113));
  bNoise(114,342) = sqrt(0.5 + 1.25/(0.125 + xConcCells(115)*xConcCells(115)*xConcCells(115)));
  bNoise(114,344) = noiseTermDAvgU;
  bNoise(114,345) = -1.*sqrt(Input::globalParameter(0)*xConcCells(114));
  bNoise(114,1043) = -1.*sqrt(xConcCells(114));
  bNoise(115,343) = sqrt(0.5 + 1.25/(0.125 + xConcCells(114)*xConcCells(114)*xConcCells(114)));
  bNoise(115,346) = noiseTermDAvgV;
  bNoise(115,347) = -1.*sqrt(Input::globalParameter(0)*xConcCells(115));
  bNoise(115,1044) = -1.*sqrt(xConcCells(115));
  bNoise(116,348) = sqrt(0.5 + 1.25/(0.125 + xConcCells(117)*xConcCells(117)*xConcCells(117)));
  bNoise(116,350) = noiseTermDAvgU;
  bNoise(116,351) = -1.*sqrt(Input::globalParameter(0)*xConcCells(116));
  bNoise(116,1049) = -1.*sqrt(xConcCells(116));
  bNoise(117,349) = sqrt(0.5 + 1.25/(0.125 + xConcCells(116)*xConcCells(116)*xConcCells(116)));
  bNoise(117,352) = noiseTermDAvgV;
  bNoise(117,353) = -1.*sqrt(Input::globalParameter(0)*xConcCells(117));
  bNoise(117,1050) = -1.*sqrt(xConcCells(117));
  bNoise(118,354) = sqrt(0.5 + 1.25/(0.125 + xConcCells(119)*xConcCells(119)*xConcCells(119)));
  bNoise(118,356) = noiseTermDAvgU;
  bNoise(118,357) = -1.*sqrt(Input::globalParameter(0)*xConcCells(118));
  bNoise(118,1055) = -1.*sqrt(xConcCells(118));
  bNoise(119,355) = sqrt(0.5 + 1.25/(0.125 + xConcCells(118)*xConcCells(118)*xConcCells(118)));
  bNoise(119,358) = noiseTermDAvgV;
  bNoise(119,359) = -1.*sqrt(Input::globalParameter(0)*xConcCells(119));
  bNoise(119,1056) = -1.*sqrt(xConcCells(119));
  bNoise(120,360) = sqrt(0.5 + 1.25/(0.125 + xConcCells(121)*xConcCells(121)*xConcCells(121)));
  bNoise(120,362) = noiseTermDAvgU;
  bNoise(120,363) = -1.*sqrt(Input::globalParameter(0)*xConcCells(120));
  bNoise(120,1061) = -1.*sqrt(xConcCells(120));
  bNoise(121,361) = sqrt(0.5 + 1.25/(0.125 + xConcCells(120)*xConcCells(120)*xConcCells(120)));
  bNoise(121,364) = noiseTermDAvgV;
  bNoise(121,365) = -1.*sqrt(Input::globalParameter(0)*xConcCells(121));
  bNoise(121,1062) = -1.*sqrt(xConcCells(121));
  bNoise(122,366) = sqrt(0.5 + 1.25/(0.125 + xConcCells(123)*xConcCells(123)*xConcCells(123)));
  bNoise(122,368) = noiseTermDAvgU;
  bNoise(122,369) = -1.*sqrt(Input::globalParameter(0)*xConcCells(122));
  bNoise(122,1067) = -1.*sqrt(xConcCells(122));
  bNoise(123,367) = sqrt(0.5 + 1.25/(0.125 + xConcCells(122)*xConcCells(122)*xConcCells(122)));
  bNoise(123,370) = noiseTermDAvgV;
  bNoise(123,371) = -1.*sqrt(Input::globalParameter(0)*xConcCells(123));
  bNoise(123,1068) = -1.*sqrt(xConcCells(123));
  bNoise(124,372) = sqrt(0.5 + 1.25/(0.125 + xConcCells(125)*xConcCells(125)*xConcCells(125)));
  bNoise(124,374) = noiseTermDAvgU;
  bNoise(124,375) = -1.*sqrt(Input::globalParameter(0)*xConcCells(124));
  bNoise(124,1073) = -1.*sqrt(xConcCells(124));
  bNoise(125,373) = sqrt(0.5 + 1.25/(0.125 + xConcCells(124)*xConcCells(124)*xConcCells(124)));
  bNoise(125,376) = noiseTermDAvgV;
  bNoise(125,377) = -1.*sqrt(Input::globalParameter(0)*xConcCells(125));
  bNoise(125,1074) = -1.*sqrt(xConcCells(125));
  bNoise(126,378) = sqrt(0.5 + 1.25/(0.125 + xConcCells(127)*xConcCells(127)*xConcCells(127)));
  bNoise(126,380) = noiseTermDAvgU;
  bNoise(126,381) = -1.*sqrt(Input::globalParameter(0)*xConcCells(126));
  bNoise(126,1079) = -1.*sqrt(xConcCells(126));
  bNoise(127,379) = sqrt(0.5 + 1.25/(0.125 + xConcCells(126)*xConcCells(126)*xConcCells(126)));
  bNoise(127,382) = noiseTermDAvgV;
  bNoise(127,383) = -1.*sqrt(Input::globalParameter(0)*xConcCells(127));
  bNoise(127,1080) = -1.*sqrt(xConcCells(127));
  bNoise(128,384) = sqrt(0.5 + 1.25/(0.125 + xConcCells(129)*xConcCells(129)*xConcCells(129)));
  bNoise(128,386) = noiseTermDAvgU;
  bNoise(128,387) = -1.*sqrt(Input::globalParameter(0)*xConcCells(128));
  bNoise(128,1085) = -1.*sqrt(xConcCells(128));
  bNoise(129,385) = sqrt(0.5 + 1.25/(0.125 + xConcCells(128)*xConcCells(128)*xConcCells(128)));
  bNoise(129,388) = noiseTermDAvgV;
  bNoise(129,389) = -1.*sqrt(Input::globalParameter(0)*xConcCells(129));
  bNoise(129,1086) = -1.*sqrt(xConcCells(129));
  bNoise(130,390) = sqrt(0.5 + 1.25/(0.125 + xConcCells(131)*xConcCells(131)*xConcCells(131)));
  bNoise(130,392) = noiseTermDAvgU;
  bNoise(130,393) = -1.*sqrt(Input::globalParameter(0)*xConcCells(130));
  bNoise(130,1091) = -1.*sqrt(xConcCells(130));
  bNoise(131,391) = sqrt(0.5 + 1.25/(0.125 + xConcCells(130)*xConcCells(130)*xConcCells(130)));
  bNoise(131,394) = noiseTermDAvgV;
  bNoise(131,395) = -1.*sqrt(Input::globalParameter(0)*xConcCells(131));
  bNoise(131,1092) = -1.*sqrt(xConcCells(131));
  bNoise(132,396) = sqrt(0.5 + 1.25/(0.125 + xConcCells(133)*xConcCells(133)*xConcCells(133)));
  bNoise(132,398) = noiseTermDAvgU;
  bNoise(132,399) = -1.*sqrt(Input::globalParameter(0)*xConcCells(132));
  bNoise(132,1097) = -1.*sqrt(xConcCells(132));
  bNoise(133,397) = sqrt(0.5 + 1.25/(0.125 + xConcCells(132)*xConcCells(132)*xConcCells(132)));
  bNoise(133,400) = noiseTermDAvgV;
  bNoise(133,401) = -1.*sqrt(Input::globalParameter(0)*xConcCells(133));
  bNoise(133,1098) = -1.*sqrt(xConcCells(133));
  bNoise(134,402) = sqrt(0.5 + 1.25/(0.125 + xConcCells(135)*xConcCells(135)*xConcCells(135)));
  bNoise(134,404) = noiseTermDAvgU;
  bNoise(134,405) = -1.*sqrt(Input::globalParameter(0)*xConcCells(134));
  bNoise(134,1103) = -1.*sqrt(xConcCells(134));
  bNoise(135,403) = sqrt(0.5 + 1.25/(0.125 + xConcCells(134)*xConcCells(134)*xConcCells(134)));
  bNoise(135,406) = noiseTermDAvgV;
  bNoise(135,407) = -1.*sqrt(Input::globalParameter(0)*xConcCells(135));
  bNoise(135,1104) = -1.*sqrt(xConcCells(135));
  bNoise(136,408) = sqrt(0.5 + 1.25/(0.125 + xConcCells(137)*xConcCells(137)*xConcCells(137)));
  bNoise(136,410) = noiseTermDAvgU;
  bNoise(136,411) = -1.*sqrt(Input::globalParameter(0)*xConcCells(136));
  bNoise(136,1109) = -1.*sqrt(xConcCells(136));
  bNoise(137,409) = sqrt(0.5 + 1.25/(0.125 + xConcCells(136)*xConcCells(136)*xConcCells(136)));
  bNoise(137,412) = noiseTermDAvgV;
  bNoise(137,413) = -1.*sqrt(Input::globalParameter(0)*xConcCells(137));
  bNoise(137,1110) = -1.*sqrt(xConcCells(137));
  bNoise(138,414) = sqrt(0.5 + 1.25/(0.125 + xConcCells(139)*xConcCells(139)*xConcCells(139)));
  bNoise(138,416) = noiseTermDAvgU;
  bNoise(138,417) = -1.*sqrt(Input::globalParameter(0)*xConcCells(138));
  bNoise(138,1115) = -1.*sqrt(xConcCells(138));
  bNoise(139,415) = sqrt(0.5 + 1.25/(0.125 + xConcCells(138)*xConcCells(138)*xConcCells(138)));
  bNoise(139,418) = noiseTermDAvgV;
  bNoise(139,419) = -1.*sqrt(Input::globalParameter(0)*xConcCells(139));
  bNoise(139,1116) = -1.*sqrt(xConcCells(139));
  bNoise(140,420) = sqrt(0.5 + 1.25/(0.125 + xConcCells(141)*xConcCells(141)*xConcCells(141)));
  bNoise(140,422) = noiseTermDAvgU;
  bNoise(140,423) = -1.*sqrt(Input::globalParameter(0)*xConcCells(140));
  bNoise(140,1121) = -1.*sqrt(xConcCells(140));
  bNoise(141,421) = sqrt(0.5 + 1.25/(0.125 + xConcCells(140)*xConcCells(140)*xConcCells(140)));
  bNoise(141,424) = noiseTermDAvgV;
  bNoise(141,425) = -1.*sqrt(Input::globalParameter(0)*xConcCells(141));
  bNoise(141,1122) = -1.*sqrt(xConcCells(141));
  bNoise(142,426) = sqrt(0.5 + 1.25/(0.125 + xConcCells(143)*xConcCells(143)*xConcCells(143)));
  bNoise(142,428) = noiseTermDAvgU;
  bNoise(142,429) = -1.*sqrt(Input::globalParameter(0)*xConcCells(142));
  bNoise(142,1127) = -1.*sqrt(xConcCells(142));
  bNoise(143,427) = sqrt(0.5 + 1.25/(0.125 + xConcCells(142)*xConcCells(142)*xConcCells(142)));
  bNoise(143,430) = noiseTermDAvgV;
  bNoise(143,431) = -1.*sqrt(Input::globalParameter(0)*xConcCells(143));
  bNoise(143,1128) = -1.*sqrt(xConcCells(143));
  bNoise(144,432) = sqrt(0.5 + 1.25/(0.125 + xConcCells(145)*xConcCells(145)*xConcCells(145)));
  bNoise(144,434) = noiseTermDAvgU;
  bNoise(144,435) = -1.*sqrt(Input::globalParameter(0)*xConcCells(144));
  bNoise(144,1133) = -1.*sqrt(xConcCells(144));
  bNoise(145,433) = sqrt(0.5 + 1.25/(0.125 + xConcCells(144)*xConcCells(144)*xConcCells(144)));
  bNoise(145,436) = noiseTermDAvgV;
  bNoise(145,437) = -1.*sqrt(Input::globalParameter(0)*xConcCells(145));
  bNoise(145,1134) = -1.*sqrt(xConcCells(145));
  bNoise(146,438) = sqrt(0.5 + 1.25/(0.125 + xConcCells(147)*xConcCells(147)*xConcCells(147)));
  bNoise(146,440) = noiseTermDAvgU;
  bNoise(146,441) = -1.*sqrt(Input::globalParameter(0)*xConcCells(146));
  bNoise(146,1139) = -1.*sqrt(xConcCells(146));
  bNoise(147,439) = sqrt(0.5 + 1.25/(0.125 + xConcCells(146)*xConcCells(146)*xConcCells(146)));
  bNoise(147,442) = noiseTermDAvgV;
  bNoise(147,443) = -1.*sqrt(Input::globalParameter(0)*xConcCells(147));
  bNoise(147,1140) = -1.*sqrt(xConcCells(147));
  bNoise(148,444) = sqrt(0.5 + 1.25/(0.125 + xConcCells(149)*xConcCells(149)*xConcCells(149)));
  bNoise(148,446) = noiseTermDAvgU;
  bNoise(148,447) = -1.*sqrt(Input::globalParameter(0)*xConcCells(148));
  bNoise(148,1145) = -1.*sqrt(xConcCells(148));
  bNoise(149,445) = sqrt(0.5 + 1.25/(0.125 + xConcCells(148)*xConcCells(148)*xConcCells(148)));
  bNoise(149,448) = noiseTermDAvgV;
  bNoise(149,449) = -1.*sqrt(Input::globalParameter(0)*xConcCells(149));
  bNoise(149,1146) = -1.*sqrt(xConcCells(149));
  bNoise(150,450) = sqrt(0.5 + 1.25/(0.125 + xConcCells(151)*xConcCells(151)*xConcCells(151)));
  bNoise(150,452) = noiseTermDAvgU;
  bNoise(150,453) = -1.*sqrt(Input::globalParameter(0)*xConcCells(150));
  bNoise(150,1151) = -1.*sqrt(xConcCells(150));
  bNoise(151,451) = sqrt(0.5 + 1.25/(0.125 + xConcCells(150)*xConcCells(150)*xConcCells(150)));
  bNoise(151,454) = noiseTermDAvgV;
  bNoise(151,455) = -1.*sqrt(Input::globalParameter(0)*xConcCells(151));
  bNoise(151,1152) = -1.*sqrt(xConcCells(151));
  bNoise(152,456) = sqrt(0.5 + 1.25/(0.125 + xConcCells(153)*xConcCells(153)*xConcCells(153)));
  bNoise(152,458) = noiseTermDAvgU;
  bNoise(152,459) = -1.*sqrt(Input::globalParameter(0)*xConcCells(152));
  bNoise(152,1157) = -1.*sqrt(xConcCells(152));
  bNoise(153,457) = sqrt(0.5 + 1.25/(0.125 + xConcCells(152)*xConcCells(152)*xConcCells(152)));
  bNoise(153,460) = noiseTermDAvgV;
  bNoise(153,461) = -1.*sqrt(Input::globalParameter(0)*xConcCells(153));
  bNoise(153,1158) = -1.*sqrt(xConcCells(153));
  bNoise(154,462) = sqrt(0.5 + 1.25/(0.125 + xConcCells(155)*xConcCells(155)*xConcCells(155)));
  bNoise(154,464) = noiseTermDAvgU;
  bNoise(154,465) = -1.*sqrt(Input::globalParameter(0)*xConcCells(154));
  bNoise(154,1163) = -1.*sqrt(xConcCells(154));
  bNoise(155,463) = sqrt(0.5 + 1.25/(0.125 + xConcCells(154)*xConcCells(154)*xConcCells(154)));
  bNoise(155,466) = noiseTermDAvgV;
  bNoise(155,467) = -1.*sqrt(Input::globalParameter(0)*xConcCells(155));
  bNoise(155,1164) = -1.*sqrt(xConcCells(155));
  bNoise(156,468) = sqrt(0.5 + 1.25/(0.125 + xConcCells(157)*xConcCells(157)*xConcCells(157)));
  bNoise(156,470) = noiseTermDAvgU;
  bNoise(156,471) = -1.*sqrt(Input::globalParameter(0)*xConcCells(156));
  bNoise(156,1169) = -1.*sqrt(xConcCells(156));
  bNoise(157,469) = sqrt(0.5 + 1.25/(0.125 + xConcCells(156)*xConcCells(156)*xConcCells(156)));
  bNoise(157,472) = noiseTermDAvgV;
  bNoise(157,473) = -1.*sqrt(Input::globalParameter(0)*xConcCells(157));
  bNoise(157,1170) = -1.*sqrt(xConcCells(157));
  bNoise(158,474) = sqrt(0.5 + 1.25/(0.125 + xConcCells(159)*xConcCells(159)*xConcCells(159)));
  bNoise(158,476) = noiseTermDAvgU;
  bNoise(158,477) = -1.*sqrt(Input::globalParameter(0)*xConcCells(158));
  bNoise(158,1175) = -1.*sqrt(xConcCells(158));
  bNoise(159,475) = sqrt(0.5 + 1.25/(0.125 + xConcCells(158)*xConcCells(158)*xConcCells(158)));
  bNoise(159,478) = noiseTermDAvgV;
  bNoise(159,479) = -1.*sqrt(Input::globalParameter(0)*xConcCells(159));
  bNoise(159,1176) = -1.*sqrt(xConcCells(159));
  bNoise(160,480) = sqrt(0.5 + 1.25/(0.125 + xConcCells(161)*xConcCells(161)*xConcCells(161)));
  bNoise(160,482) = noiseTermDAvgU;
  bNoise(160,483) = -1.*sqrt(Input::globalParameter(0)*xConcCells(160));
  bNoise(160,1181) = -1.*sqrt(xConcCells(160));
  bNoise(161,481) = sqrt(0.5 + 1.25/(0.125 + xConcCells(160)*xConcCells(160)*xConcCells(160)));
  bNoise(161,484) = noiseTermDAvgV;
  bNoise(161,485) = -1.*sqrt(Input::globalParameter(0)*xConcCells(161));
  bNoise(161,1182) = -1.*sqrt(xConcCells(161));
  bNoise(162,486) = sqrt(0.5 + 1.25/(0.125 + xConcCells(163)*xConcCells(163)*xConcCells(163)));
  bNoise(162,488) = noiseTermDAvgU;
  bNoise(162,489) = -1.*sqrt(Input::globalParameter(0)*xConcCells(162));
  bNoise(162,1187) = -1.*sqrt(xConcCells(162));
  bNoise(163,487) = sqrt(0.5 + 1.25/(0.125 + xConcCells(162)*xConcCells(162)*xConcCells(162)));
  bNoise(163,490) = noiseTermDAvgV;
  bNoise(163,491) = -1.*sqrt(Input::globalParameter(0)*xConcCells(163));
  bNoise(163,1188) = -1.*sqrt(xConcCells(163));
  bNoise(164,492) = sqrt(0.5 + 1.25/(0.125 + xConcCells(165)*xConcCells(165)*xConcCells(165)));
  bNoise(164,494) = noiseTermDAvgU;
  bNoise(164,495) = -1.*sqrt(Input::globalParameter(0)*xConcCells(164));
  bNoise(164,1193) = -1.*sqrt(xConcCells(164));
  bNoise(165,493) = sqrt(0.5 + 1.25/(0.125 + xConcCells(164)*xConcCells(164)*xConcCells(164)));
  bNoise(165,496) = noiseTermDAvgV;
  bNoise(165,497) = -1.*sqrt(Input::globalParameter(0)*xConcCells(165));
  bNoise(165,1194) = -1.*sqrt(xConcCells(165));
  bNoise(166,498) = sqrt(0.5 + 1.25/(0.125 + xConcCells(167)*xConcCells(167)*xConcCells(167)));
  bNoise(166,500) = noiseTermDAvgU;
  bNoise(166,501) = -1.*sqrt(Input::globalParameter(0)*xConcCells(166));
  bNoise(166,1199) = -1.*sqrt(xConcCells(166));
  bNoise(167,499) = sqrt(0.5 + 1.25/(0.125 + xConcCells(166)*xConcCells(166)*xConcCells(166)));
  bNoise(167,502) = noiseTermDAvgV;
  bNoise(167,503) = -1.*sqrt(Input::globalParameter(0)*xConcCells(167));
  bNoise(167,1200) = -1.*sqrt(xConcCells(167));
  bNoise(168,504) = sqrt(0.5 + 1.25/(0.125 + xConcCells(169)*xConcCells(169)*xConcCells(169)));
  bNoise(168,506) = noiseTermDAvgU;
  bNoise(168,507) = -1.*sqrt(Input::globalParameter(0)*xConcCells(168));
  bNoise(168,1205) = -1.*sqrt(xConcCells(168));
  bNoise(169,505) = sqrt(0.5 + 1.25/(0.125 + xConcCells(168)*xConcCells(168)*xConcCells(168)));
  bNoise(169,508) = noiseTermDAvgV;
  bNoise(169,509) = -1.*sqrt(Input::globalParameter(0)*xConcCells(169));
  bNoise(169,1206) = -1.*sqrt(xConcCells(169));
  bNoise(170,510) = sqrt(0.5 + 1.25/(0.125 + xConcCells(171)*xConcCells(171)*xConcCells(171)));
  bNoise(170,512) = noiseTermDAvgU;
  bNoise(170,513) = -1.*sqrt(Input::globalParameter(0)*xConcCells(170));
  bNoise(170,1211) = -1.*sqrt(xConcCells(170));
  bNoise(171,511) = sqrt(0.5 + 1.25/(0.125 + xConcCells(170)*xConcCells(170)*xConcCells(170)));
  bNoise(171,514) = noiseTermDAvgV;
  bNoise(171,515) = -1.*sqrt(Input::globalParameter(0)*xConcCells(171));
  bNoise(171,1212) = -1.*sqrt(xConcCells(171));
  bNoise(172,516) = sqrt(0.5 + 1.25/(0.125 + xConcCells(173)*xConcCells(173)*xConcCells(173)));
  bNoise(172,518) = noiseTermDAvgU;
  bNoise(172,519) = -1.*sqrt(Input::globalParameter(0)*xConcCells(172));
  bNoise(172,1217) = -1.*sqrt(xConcCells(172));
  bNoise(173,517) = sqrt(0.5 + 1.25/(0.125 + xConcCells(172)*xConcCells(172)*xConcCells(172)));
  bNoise(173,520) = noiseTermDAvgV;
  bNoise(173,521) = -1.*sqrt(Input::globalParameter(0)*xConcCells(173));
  bNoise(173,1218) = -1.*sqrt(xConcCells(173));
  bNoise(174,522) = sqrt(0.5 + 1.25/(0.125 + xConcCells(175)*xConcCells(175)*xConcCells(175)));
  bNoise(174,524) = noiseTermDAvgU;
  bNoise(174,525) = -1.*sqrt(Input::globalParameter(0)*xConcCells(174));
  bNoise(174,1223) = -1.*sqrt(xConcCells(174));
  bNoise(175,523) = sqrt(0.5 + 1.25/(0.125 + xConcCells(174)*xConcCells(174)*xConcCells(174)));
  bNoise(175,526) = noiseTermDAvgV;
  bNoise(175,527) = -1.*sqrt(Input::globalParameter(0)*xConcCells(175));
  bNoise(175,1224) = -1.*sqrt(xConcCells(175));
  bNoise(176,528) = sqrt(0.5 + 1.25/(0.125 + xConcCells(177)*xConcCells(177)*xConcCells(177)));
  bNoise(176,530) = noiseTermDAvgU;
  bNoise(176,531) = -1.*sqrt(Input::globalParameter(0)*xConcCells(176));
  bNoise(176,1229) = -1.*sqrt(xConcCells(176));
  bNoise(177,529) = sqrt(0.5 + 1.25/(0.125 + xConcCells(176)*xConcCells(176)*xConcCells(176)));
  bNoise(177,532) = noiseTermDAvgV;
  bNoise(177,533) = -1.*sqrt(Input::globalParameter(0)*xConcCells(177));
  bNoise(177,1230) = -1.*sqrt(xConcCells(177));
  bNoise(178,534) = sqrt(0.5 + 1.25/(0.125 + xConcCells(179)*xConcCells(179)*xConcCells(179)));
  bNoise(178,536) = noiseTermDAvgU;
  bNoise(178,537) = -1.*sqrt(Input::globalParameter(0)*xConcCells(178));
  bNoise(178,1235) = -1.*sqrt(xConcCells(178));
  bNoise(179,535) = sqrt(0.5 + 1.25/(0.125 + xConcCells(178)*xConcCells(178)*xConcCells(178)));
  bNoise(179,538) = noiseTermDAvgV;
  bNoise(179,539) = -1.*sqrt(Input::globalParameter(0)*xConcCells(179));
  bNoise(179,1236) = -1.*sqrt(xConcCells(179));
  bNoise(180,540) = sqrt(0.5 + 1.25/(0.125 + xConcCells(181)*xConcCells(181)*xConcCells(181)));
  bNoise(180,542) = noiseTermDAvgU;
  bNoise(180,543) = -1.*sqrt(Input::globalParameter(0)*xConcCells(180));
  bNoise(180,1241) = -1.*sqrt(xConcCells(180));
  bNoise(181,541) = sqrt(0.5 + 1.25/(0.125 + xConcCells(180)*xConcCells(180)*xConcCells(180)));
  bNoise(181,544) = noiseTermDAvgV;
  bNoise(181,545) = -1.*sqrt(Input::globalParameter(0)*xConcCells(181));
  bNoise(181,1242) = -1.*sqrt(xConcCells(181));
  bNoise(182,546) = sqrt(0.5 + 1.25/(0.125 + xConcCells(183)*xConcCells(183)*xConcCells(183)));
  bNoise(182,548) = noiseTermDAvgU;
  bNoise(182,549) = -1.*sqrt(Input::globalParameter(0)*xConcCells(182));
  bNoise(182,1247) = -1.*sqrt(xConcCells(182));
  bNoise(183,547) = sqrt(0.5 + 1.25/(0.125 + xConcCells(182)*xConcCells(182)*xConcCells(182)));
  bNoise(183,550) = noiseTermDAvgV;
  bNoise(183,551) = -1.*sqrt(Input::globalParameter(0)*xConcCells(183));
  bNoise(183,1248) = -1.*sqrt(xConcCells(183));
  bNoise(184,552) = sqrt(0.5 + 1.25/(0.125 + xConcCells(185)*xConcCells(185)*xConcCells(185)));
  bNoise(184,554) = noiseTermDAvgU;
  bNoise(184,555) = -1.*sqrt(Input::globalParameter(0)*xConcCells(184));
  bNoise(184,1253) = -1.*sqrt(xConcCells(184));
  bNoise(185,553) = sqrt(0.5 + 1.25/(0.125 + xConcCells(184)*xConcCells(184)*xConcCells(184)));
  bNoise(185,556) = noiseTermDAvgV;
  bNoise(185,557) = -1.*sqrt(Input::globalParameter(0)*xConcCells(185));
  bNoise(185,1254) = -1.*sqrt(xConcCells(185));
  bNoise(186,558) = sqrt(0.5 + 1.25/(0.125 + xConcCells(187)*xConcCells(187)*xConcCells(187)));
  bNoise(186,560) = noiseTermDAvgU;
  bNoise(186,561) = -1.*sqrt(Input::globalParameter(0)*xConcCells(186));
  bNoise(186,1259) = -1.*sqrt(xConcCells(186));
  bNoise(187,559) = sqrt(0.5 + 1.25/(0.125 + xConcCells(186)*xConcCells(186)*xConcCells(186)));
  bNoise(187,562) = noiseTermDAvgV;
  bNoise(187,563) = -1.*sqrt(Input::globalParameter(0)*xConcCells(187));
  bNoise(187,1260) = -1.*sqrt(xConcCells(187));
  bNoise(188,564) = sqrt(0.5 + 1.25/(0.125 + xConcCells(189)*xConcCells(189)*xConcCells(189)));
  bNoise(188,566) = noiseTermDAvgU;
  bNoise(188,567) = -1.*sqrt(Input::globalParameter(0)*xConcCells(188));
  bNoise(188,1265) = -1.*sqrt(xConcCells(188));
  bNoise(189,565) = sqrt(0.5 + 1.25/(0.125 + xConcCells(188)*xConcCells(188)*xConcCells(188)));
  bNoise(189,568) = noiseTermDAvgV;
  bNoise(189,569) = -1.*sqrt(Input::globalParameter(0)*xConcCells(189));
  bNoise(189,1266) = -1.*sqrt(xConcCells(189));
  bNoise(190,570) = sqrt(0.5 + 1.25/(0.125 + xConcCells(191)*xConcCells(191)*xConcCells(191)));
  bNoise(190,572) = noiseTermDAvgU;
  bNoise(190,573) = -1.*sqrt(Input::globalParameter(0)*xConcCells(190));
  bNoise(190,1271) = -1.*sqrt(xConcCells(190));
  bNoise(191,571) = sqrt(0.5 + 1.25/(0.125 + xConcCells(190)*xConcCells(190)*xConcCells(190)));
  bNoise(191,574) = noiseTermDAvgV;
  bNoise(191,575) = -1.*sqrt(Input::globalParameter(0)*xConcCells(191));
  bNoise(191,1272) = -1.*sqrt(xConcCells(191));
  bNoise(192,576) = sqrt(0.5 + 1.25/(0.125 + xConcCells(193)*xConcCells(193)*xConcCells(193)));
  bNoise(192,578) = noiseTermDAvgU;
  bNoise(192,579) = -1.*sqrt(Input::globalParameter(0)*xConcCells(192));
  bNoise(192,1277) = -1.*sqrt(xConcCells(192));
  bNoise(193,577) = sqrt(0.5 + 1.25/(0.125 + xConcCells(192)*xConcCells(192)*xConcCells(192)));
  bNoise(193,580) = noiseTermDAvgV;
  bNoise(193,581) = -1.*sqrt(Input::globalParameter(0)*xConcCells(193));
  bNoise(193,1278) = -1.*sqrt(xConcCells(193));
  bNoise(194,582) = sqrt(0.5 + 1.25/(0.125 + xConcCells(195)*xConcCells(195)*xConcCells(195)));
  bNoise(194,584) = noiseTermDAvgU;
  bNoise(194,585) = -1.*sqrt(Input::globalParameter(0)*xConcCells(194));
  bNoise(194,1283) = -1.*sqrt(xConcCells(194));
  bNoise(195,583) = sqrt(0.5 + 1.25/(0.125 + xConcCells(194)*xConcCells(194)*xConcCells(194)));
  bNoise(195,586) = noiseTermDAvgV;
  bNoise(195,587) = -1.*sqrt(Input::globalParameter(0)*xConcCells(195));
  bNoise(195,1284) = -1.*sqrt(xConcCells(195));
  bNoise(196,588) = sqrt(0.5 + 1.25/(0.125 + xConcCells(197)*xConcCells(197)*xConcCells(197)));
  bNoise(196,590) = noiseTermDAvgU;
  bNoise(196,591) = -1.*sqrt(Input::globalParameter(0)*xConcCells(196));
  bNoise(196,1289) = -1.*sqrt(xConcCells(196));
  bNoise(197,589) = sqrt(0.5 + 1.25/(0.125 + xConcCells(196)*xConcCells(196)*xConcCells(196)));
  bNoise(197,592) = noiseTermDAvgV;
  bNoise(197,593) = -1.*sqrt(Input::globalParameter(0)*xConcCells(197));
  bNoise(197,1290) = -1.*sqrt(xConcCells(197));
  bNoise(198,594) = sqrt(0.5 + 1.25/(0.125 + xConcCells(199)*xConcCells(199)*xConcCells(199)));
  bNoise(198,596) = noiseTermDAvgU;
  bNoise(198,597) = -1.*sqrt(Input::globalParameter(0)*xConcCells(198));
  bNoise(198,1295) = -1.*sqrt(xConcCells(198));
  bNoise(199,595) = sqrt(0.5 + 1.25/(0.125 + xConcCells(198)*xConcCells(198)*xConcCells(198)));
  bNoise(199,598) = noiseTermDAvgV;
  bNoise(199,599) = -1.*sqrt(Input::globalParameter(0)*xConcCells(199));
  bNoise(199,1296) = -1.*sqrt(xConcCells(199));
}

void ChemicalLangevinComputeIncrementProfiling::computeMatrixVectorProduct
(
  Array<double,2>& bNoise,
  Array<double,1>& randomGaussianVector,
  Array<double,1>& noiseIncrementVector
)
{
  // Compute the matrix-vector product noiseIncrementVector = BNoise * randomGaussianVector.
//  noiseIncrementVector.resize(nStot);

  // Matrix-vector multiplication using Blitz++ library
//  blitz::firstIndex ii;
//  blitz::secondIndex jj;
//  noiseIncrementVector = blitz::sum(bNoise(ii,jj) * randomGaussianVector(jj), jj);

  // BLAS matrix-vector multiplication
  // Note: We use the optimized BLAS library from ATLAS
  // (Automatically Tuned Linear Algebra Software) which is
  // architecture-specific.
  // See http://math-atlas.sourceforge.net/

  const int m = bNoise.extent(blitz::firstDim);
  const int n = bNoise.extent(blitz::secondDim);
  const double alpha = 1.0;
  const double beta = 0.0;
  const double *bNoisePointer = bNoise.data();
  const double *randomGaussianVectorPointer = randomGaussianVector.data();
  double *noiseIncrementVectorPointer = noiseIncrementVector.data();
  const int lda = n;
  const int incX = 1;
  const int incY = 1;
  cblas_dgemv(CblasRowMajor,CblasNoTrans,m,n,alpha,bNoisePointer,lda,
              randomGaussianVectorPointer,incX,beta,
              noiseIncrementVectorPointer,incY);
  //cout << noiseIncrementVector << endl;
}

void ChemicalLangevinComputeIncrementProfiling::addNoiseIncrement
(
  Array<double,1>& xConcCells,
  Array<double,1>& xConcMilieu,
  const int nCells,
  const double volumeCell,
  const double volumeExt,
  const double chemicalLangevinTimeStep
)
{
  // Add the noise increments to the xConc vectors.
  double sqrtTimeStep = sqrt(chemicalLangevinTimeStep);
  xConcCells += (sqrt(nMConversion/volumeCell))*noiseIncrementVector(blitz::Range(0,nSCellstot-1))*sqrtTimeStep;
  xConcMilieu += (sqrt(nMConversion/volumeExt))*noiseIncrementVector(blitz::Range(nSCellstot,nStot-1))*sqrtTimeStep;
}


void ChemicalLangevinComputeIncrementProfiling::chemicalLangevinComputeIncrement
(
  Array<double,1>& xConcCells,
  Array<double,1>& xConcMilieu,
  const int nCells,
  const double volumeCell,
  const double volumeExt,
  const double chemicalLangevinTimeStep
)
{

  nSCellstot = nCells*nSCell;
  nRtot = 2.*(nCells*(nRCell+nRCellMilieu) + nRMilieu);
  nStot = nCells*nSCell + nSMilieu;

  // Computing the population average of U and V
  Array<double,1> dummy1( xConcCells( blitz::Range( 0, blitz::toEnd, nSCell ) ) );
  Array<double,1> dummy2( xConcCells( blitz::Range( 1, blitz::toEnd, nSCell ) ) );
  double avgU = computeMean( dummy1 );
  double avgV = computeMean( dummy2 );

  computeDeterministicIncrement
  (
    xConcCells,
    xConcMilieu,
    nCells,
    volumeCell,
    volumeExt,
    chemicalLangevinTimeStep,
    avgU,
    avgV
  );


  Array<double,1> randomGaussianVector(nRtot);
  generateRandomVector(randomGaussianVector);

  computeBNoise(bNoise,
                xConcCells,
                xConcMilieu,
                avgU,
                avgV
                );

  computeMatrixVectorProduct
  (
    bNoise,
    randomGaussianVector,
    noiseIncrementVector
  );

  addNoiseIncrement
  (
    xConcCells,
    xConcMilieu,
    nCells,
    volumeCell,
    volumeExt,
    chemicalLangevinTimeStep
  );

  // Correct the xConc for negative values
  for (int i=0; i<nSCellstot; ++i)
  {
    xConcCells(i) = abs(xConcCells(i));
  }
  for (int i=0; i<nSMilieu; ++i)
  {
    xConcMilieu(i) = abs(xConcMilieu(i));
  }

}

#endif // USE_CHEMICAL_LANGEVIN
