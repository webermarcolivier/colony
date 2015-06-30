/***************************************************************************//**
 * Project: Colony
 *
 * \file    CellMilieuChemicalSystem.cpp
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

#include "CellMilieuChemicalSystem.h"

//------------------------------------------------------------------------------

CellMilieuChemicalSystem::CellMilieuChemicalSystem()
  : milieuPtr_(0),
    cellPtr_(0),
    nTotalSpecies_(0),
    nSpeciesCell_(0),
    nSpeciesMilieu_(0),
    nChannels_(0),
    stoichMatrix_(), //Rem: the Array is 0-initialized, its data cannot be accessed
    computePropensitiesInterCellularPtr_(0),
    computeTimeDependentPropensitiesInterCellularPtr_(0),
    propensities_() //Rem: the Array is 0-initialized, its data cannot be accessed
{}

//------------------------------------------------------------------------------

CellMilieuChemicalSystem::~CellMilieuChemicalSystem()
{
  stoichMatrix_.free();
  propensities_.free();
}

//------------------------------------------------------------------------------

void CellMilieuChemicalSystem::initialize
(
  const CellMilieuChemicalSystemInitParam& p
)
{
  milieuPtr_ = p.milieuPtr_;
  cellPtr_   = p.cellPtr_;

  nSpeciesMilieu_ = milieuPtr_->getNSpecies();
  nSpeciesCell_   =   cellPtr_->getNSpecies();
  nTotalSpecies_  = nSpeciesMilieu_ + nSpeciesCell_;
  nChannels_      = p.nChannelsCellMilieu_;

  int nrows = p.stoichMatrixCellMilieu_.rows();
  int ncols = p.stoichMatrixCellMilieu_.cols();
  if (nrows != nChannels_ || ncols != nTotalSpecies_)
  {
    cout << "ERROR: class CellMilieuChemicalSystem initialize(),"
      " stoichMatrix array does not have the correct size:\n";
    cout << "nChannels x nSpecies = " << nChannels_ << " x "
         << nTotalSpecies_ << "\n";
    cout << "dim(stoichMatrix) = " << nrows << " x " << ncols << endl;
    exit(1);
  }
  stoichMatrix_.resize(nrows,ncols);
  stoichMatrix_ = p.stoichMatrixCellMilieu_;

  propensities_.resize(nChannels_);
  propensities_ = 0.0;

  #ifndef TIME_DEPENDENT_PROPENSITIES
    computePropensitiesInterCellularPtr_ = p.computePropensitiesInterCellularPtr_;
  #else
    computeTimeDependentPropensitiesInterCellularPtr_ =
                            p.computeTimeDependentPropensitiesInterCellularPtr_;
  #endif // TIME_DEPENDENT_PROPENSITIES
}

//------------------------------------------------------------------------------

void CellMilieuChemicalSystem::setReferencePropensities(Array<double,1> refArray)
{
  Array<double,1> arrayCopy( propensities_.copy() );
  propensities_.reference(refArray);
  propensities_ = arrayCopy;
}

//------------------------------------------------------------------------------

void CellMilieuChemicalSystem::computePropensities()
{
  computePropensitiesInterCellularPtr_
  (
    cellPtr_->getX(),
    milieuPtr_->getX(),
    cellPtr_->getVolume0(),
    milieuPtr_->getTimeDependentVolume(0.0),
    propensities_
  );
}

//------------------------------------------------------------------------------

void CellMilieuChemicalSystem::computePropensities(const double time)
{
  // UPDATE: We changed the definition of the volume of the milieu, but this is kind of
  // confusing because we have to use the getTimeDependentVolume(time) method to get
  // the value of the milieu volume, even in the case of time independent propensities.
  // TODO: find a better scheme for definition of the volumes of cells and milieu
  // that is compatible with both time dependent and time independent propensities.

  computeTimeDependentPropensitiesInterCellularPtr_
  (
    cellPtr_->getX(),
    milieuPtr_->getX(),    
    cellPtr_->getTimeDependentVolume(time),
    cellPtr_->getVolume0(),
    milieuPtr_->getTimeDependentVolume(time),
    propensities_
  );
}

//------------------------------------------------------------------------------

void CellMilieuChemicalSystem::applyReaction(int mu)
{
  cellPtr_->applyInterCellularReaction(
                stoichMatrix_( mu, Range( 0, nSpeciesCell_ - 1 ) )
                                      );
  milieuPtr_->applyInterCellularReaction(
                stoichMatrix_( mu, Range( nSpeciesCell_, nTotalSpecies_ - 1 ) )
                                        );

  #ifdef X_POSITIVITY_CHECK
  checkPositivity();
  #endif // X_POSITIVITY_CHECK
}

//------------------------------------------------------------------------------

#ifdef X_POSITIVITY_CHECK
void CellMilieuChemicalSystem::checkPositivity() const
{
  int i;
  for (i=0; i<nSpeciesCell_; i++)
  {
    if ( (cellPtr_->getX())(i) < 0)
    {
      cout << "WARNING: class CellMilieuChemicalSystem, function checkPositivity(),"
              " state vector x_ of cell contains a negative number of molecules: "
              "x(" << i << ") = " << (cellPtr_->getX())(i) << endl;
    }
  }

  for (i=0; i<nSpeciesMilieu_; i++)
  {
    if ( (milieuPtr_->getX())(i) < 0)
    {
      cout << "WARNING: class CellMilieuChemicalSystem, function checkPositivity(),"
              " state vector x_ of milieu contains a negative number of molecules: "
              "x(" << i << ") = " << (milieuPtr_->getX())(i) << endl;
    }
  }
}
#endif //X_POSITIVITY_CHECK

//------------------------------------------------------------------------------

CellMilieuChemicalSystem& CellMilieuChemicalSystem::operator =
(
  const CellMilieuChemicalSystem& b
)
{
  // Test for self-assignment
  if (&b == this)
  {
    // Do nothing and just return myself.
    return *this;
  }

  copy(b);

  return *this;
}

//------------------------------------------------------------------------------

void CellMilieuChemicalSystem::copy(const CellMilieuChemicalSystem& b)
{
  milieuPtr_ = b.milieuPtr_;

  // Warning: We don't want to make the cellMilieu object point to the same cell!!!
  //  cellPtr_   = b.cellPtr_;

  nSpeciesMilieu_ = b.nSpeciesMilieu_;
  nSpeciesCell_   = b.nSpeciesCell_;
  nTotalSpecies_  = b.nTotalSpecies_;
  nChannels_      = b.nChannels_;

  int nrows = b.stoichMatrix_.rows();
  int ncols = b.stoichMatrix_.cols();
  stoichMatrix_.resize(nrows,ncols);
  stoichMatrix_ = b.stoichMatrix_;

  propensities_.resize(b.propensities_.size());
  propensities_ = b.propensities_;

  computePropensitiesInterCellularPtr_ = b.computePropensitiesInterCellularPtr_;
  computeTimeDependentPropensitiesInterCellularPtr_ =
                            b.computeTimeDependentPropensitiesInterCellularPtr_;
}

//------------------------------------------------------------------------------









