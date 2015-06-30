/***************************************************************************//**
 * Project: Colony
 *
 * \file    Cell.cpp
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

#include "Cell.h"

//------------------------------------------------------------------------------

/*
 * The definition/initialization of the static members have to go in the .cpp
 * file. If put in the .h file, then all classes that include the .h file will
 * try to declare the static member which can be decclared only once.
 */

int Cell::nCellObjects_ = 0;

//------------------------------------------------------------------------------

Cell::Cell()
  : cellCycleDeterministic_(0.0),
    cellCycleGamma_(0.0),
    nTotalChannels_(0),
    totalPropensities_(), //Rem: the Array is 0-initialized, its data cannot be accessed
    dnaSpecies_() //Rem: the Array is 0-initialized, its data cannot be accessed
{
  nCellObjects_++;
  GraphicsCellComposite::setCellIndex(-1);
}

//------------------------------------------------------------------------------

Cell::~Cell()
{
  nCellObjects_--;
  totalPropensities_.free();
  dnaSpecies_.free();
}

//------------------------------------------------------------------------------

void Cell::initialize (CellInitParam& p)
{
  CellBase::initialize(p.cellBaseInitParam_);
  GraphicsCellComposite::initialize(p.graphicsCellParam_);

  // Initialize the parameters for the cell cycle distribution.
  cellCycleDeterministic_ = p.cellCycleDeterministic_;
  cellCycleGamma_         = p.cellCycleGamma_;

  // Initialize the cell cycle initial phase and duration.
  state_.cellCycleDuration_ =
    RandomNumberGenerator::getCellCycle(cellCycleDeterministic_,cellCycleGamma_);
  state_.timePreviousDivision_ = - p.cellCyclephase_*state_.cellCycleDuration_;

  // Compute the initial xConc
  state_.xConc_ = CellBase::getXconc(0.0);

  // Test if the array dnaSpecies have correct size.
  int nSpecies = p.cellBaseInitParam_.chemicalSystemInitParam_.nSpecies_;
  if (nSpecies != p.dnaSpecies_.rows() ||
      nSpecies != p.dnaSpecies_.cols())
  {
    cout << "ERROR: class Cell initialize(), dnaSpecies array does not have"
      " the correct size:\n";
    cout << "nSpecies = " << nSpecies << "\n";
    cout << "dimensions(dnaSpecies) = " << p.dnaSpecies_.rows()
         << " x " << p.dnaSpecies_.cols() << endl;
    exit(1);
  }

  dnaSpecies_.resize(nSpecies,nSpecies);
  dnaSpecies_ = p.dnaSpecies_;

  int nChannels = p.cellBaseInitParam_.chemicalSystemInitParam_.nChannels_;
  nTotalChannels_ = nChannels +
                    p.cellMilieuChemicalSystemInitParam_.nChannelsCellMilieu_;

  totalPropensities_.resize(nTotalChannels_);

  CellBase::setReferencePropensities(totalPropensities_(Range(0, nChannels-1)));

  // Initiate the CellMilieuChemicalSystem.
  p.cellMilieuChemicalSystemInitParam_.cellPtr_ = this;
  cellMilieuChemicalSystem_.initialize(p.cellMilieuChemicalSystemInitParam_);


  cellMilieuChemicalSystem_.setReferencePropensities
  (
    totalPropensities_( Range(nChannels, nTotalChannels_ - 1) )
  );

  totalPropensities_ = 0.0;
}

//------------------------------------------------------------------------------

void Cell::copy(const Cell& c2)
{
  CellBase::copy(c2);
  GraphicsCellComposite::copy(c2);

  cellCycleDeterministic_ = c2.cellCycleDeterministic_;
  cellCycleGamma_         = c2.cellCycleGamma_;

  int dim1 = c2.dnaSpecies_.extent(firstDim);
  int dim2 = c2.dnaSpecies_.extent(secondDim);
  dnaSpecies_.resize(dim1,dim2);
  dnaSpecies_ = c2.dnaSpecies_;

  nTotalChannels_ = c2.nTotalChannels_;

  int dim = c2.totalPropensities_.size();
  totalPropensities_.resize(dim);

  int nChannels = c2.chemicalSystem_.nChannels_;
  CellBase::setReferencePropensities(totalPropensities_(Range(0, nChannels-1)));

  // CellMilieuChemicalSystem assigment operator.
  cellMilieuChemicalSystem_ = c2.cellMilieuChemicalSystem_;
  // Make the cellMilieu object point to myself.
  cellMilieuChemicalSystem_.cellPtr_ = this;


  cellMilieuChemicalSystem_.setReferencePropensities
  (
    totalPropensities_( Range(nChannels, nTotalChannels_ - 1) )
  );

  totalPropensities_ = c2.totalPropensities_;
}

//------------------------------------------------------------------------------

Cell::Cell(const Cell& c2)
  : GraphicsCellComposite(),
    CellBase(),
    nTotalChannels_(0),
    totalPropensities_(), //Rem: the Array is 0-initialized, its data cannot be accessed
    dnaSpecies_() //Rem: the Array is 0-initialized, its data cannot be accessed
{

  nCellObjects_++;

  copy(c2);
}

//------------------------------------------------------------------------------

Cell& Cell::operator = (const Cell& c2)
{
  // Test for self-assignment
  if (&c2 == this)
  {
    // Do nothing and just return myself.
    return *this;
  }

  copy(c2);

  return *this;
}

//------------------------------------------------------------------------------

ostream& operator<< (ostream& out, const Cell& c)
{
  operator<<(out,(CellBase&)c);
  /*out << "nChannels cell-milieu = "
      << c.cellMilieuChemicalSystem_.getNChannels() << '\n';*/
  /*out << "Dimensions(dnaSpecies) = "
      << c.dnaSpecies_.extent(firstDim) << " x "
      << c.dnaSpecies_.extent(secondDim) << "\n";*/
  return out;
}

//------------------------------------------------------------------------------

void Cell::setCellIndex(int i)
{
  GraphicsCellComposite::setCellIndex(i);
}

//------------------------------------------------------------------------------

void Cell::setReferencePropensities(Array<double,1> refArray)
{
  totalPropensities_.reference(refArray);

  int nChannels = CellBase::getNChannels();
  CellBase::setReferencePropensities(totalPropensities_(Range(0, nChannels-1)));

  int nChannelsCellMilieu = cellMilieuChemicalSystem_.getNChannels();
  cellMilieuChemicalSystem_.setReferencePropensities
  (
    totalPropensities_( Range(nChannels, nChannels + nChannelsCellMilieu - 1) )
  );
}

//------------------------------------------------------------------------------

void Cell::setPosition(TinyVector<double,3> position)
{
  CellBase::setPosition(position);
  GraphicsCellComposite::setPosition(position(0), position(1), position(2));
}

//------------------------------------------------------------------------------

void Cell::setAngle(double angle)
{
  CellBase::setAngle(angle);
  GraphicsCellComposite::setAngle(angle);
}

//------------------------------------------------------------------------------

double Cell::getAngle() const
{
  return CellBase::getAngle();
}



//------------------------------------------------------------------------------

void Cell::updateGraphicsCell(const double time)
{
  TinyVector<double,3> cellDimensions = CellBase::getTimeDependentCellDimensions(time);
  GraphicsCellComposite::setCellHeight( cellDimensions(1) );
  GraphicsCellComposite::setCellLength( cellDimensions(2) );
  GraphicsCellComposite::setAngle(CellBase::state_.angle_);
  TinyVector<double,3> position(CellBase::state_.position_);
  GraphicsCellComposite::setPosition(position(0), position(1), position(2));
}

//------------------------------------------------------------------------------

void Cell::applyReaction(int mu)
{
  if ( mu < chemicalSystem_.nChannels_ )
  {
    CellBase::applyReaction(mu);

  } else if ( mu < nTotalChannels_ )
  {
    cellMilieuChemicalSystem_.applyReaction( mu - chemicalSystem_.nChannels_ );

  } else {
    cout << "ERROR: class Cell, applyReaction(): channel number mu is bigger"
    " than the total number of channels (internal + cell-milieu).\n"
    << "mu = " << mu << '\n'
    << "nTotalChannels = " << nTotalChannels_ << endl;
    exit(1);
  }
}

//------------------------------------------------------------------------------

void Cell::computePropensities()
{
  CellBase::computePropensities();
  cellMilieuChemicalSystem_.computePropensities();
}

//------------------------------------------------------------------------------

void Cell::computePropensities(double time)
{
  // This will use the time-dependent procedure of CellBase to compute the
  // propensities.

  CellBase::computePropensities(time);

  // The diffusion reactions are time-dependent and are computed by a different
  // function.
  cellMilieuChemicalSystem_.computePropensities(time);
}

//------------------------------------------------------------------------------

void Cell::computePropensitiesCell()
{
  CellBase::computePropensities();
}

//------------------------------------------------------------------------------

void Cell::computePropensitiesCell(const double time)
{
  CellBase::computePropensities(time);
}

//------------------------------------------------------------------------------

void Cell::computePropensitiesCellMilieu()
{
  cellMilieuChemicalSystem_.computePropensities();
}

//------------------------------------------------------------------------------

void Cell::computePropensitiesCellMilieu(const double time)
{
  cellMilieuChemicalSystem_.computePropensities(time);
}

//------------------------------------------------------------------------------

void Cell::detachProteinsFromDna()
{
  int i, j;
  for (i=0; i<chemicalSystem_.nSpecies_; i++)
  {
    if (dnaSpecies_(i,i) == 1)
      // See definition of dnaSpecies_ in header file of the Cell class.
    {
      // For protein-DNA complexes, we detach the proteins from the DNA.
      int nProteinDnaComplex = state_.x_(i);
      state_.x_(i) = 0;

      for (j=0; j<chemicalSystem_.nSpecies_; j++)
      {
        if (i != j)
        {
          state_.x_(j) += nProteinDnaComplex * dnaSpecies_(i,j);
            // 'dnaSpecies' defines the stoichiometric composition of the
            // protein-DNA complex.
        }
      }
    }
    // For protein and DNA molecule species, do nothing.
  }
}

//------------------------------------------------------------------------------

Cell Cell::duplicate()
{
  detachProteinsFromDna();

  // Make a copy of the mother cell.
  Cell daughterCell(*this);

  // Make two copies of the mother cell's state.
  Array<int,1> x1(state_.x_.copy());
  Array<int,1> x2(state_.x_.copy());

  // Distribute proteins to mother and daughter cells following a binomial
  // distribution.
  int j;
  for (j=0; j<chemicalSystem_.nSpecies_; j++)
  {
    if ( dnaSpecies_(j,j) == -1 )
    {
      /* This species is genetic material, copy it to both cells. The daughter
         cell is already a copy of the mother cell, do nothing. */

    } else {
      /* This species is a protein or protein-DNA complex, distribute it
         between both cells. */

      int binomial = RandomNumberGenerator::getBinomial(0.5, unsigned(x1(j)) );

      x1(j) = binomial;
      x2(j) = x2(j) - binomial;
    }
  }

  // Write the binomial distribution in the state vector of both cells.
  state_.x_              = x1;
  daughterCell.state_.x_ = x2;

  // Update the time of the last division and the new cell cycle duration
  // of the mother cell.
  state_.timePreviousDivision_ += state_.cellCycleDuration_;
  state_.cellCycleDuration_ =
    RandomNumberGenerator::getCellCycle(cellCycleDeterministic_,cellCycleGamma_);

  // Update the time of the last division and the new cell cycle duration
  // of the daughter cell.
  daughterCell.state_.timePreviousDivision_ = state_.timePreviousDivision_;
  daughterCell.state_.cellCycleDuration_ =
       RandomNumberGenerator::getCellCycle(daughterCell.cellCycleDeterministic_,
                                                  daughterCell.cellCycleGamma_);

  // Set position of the cells
  TinyVector<double,3> initialPos = getPosition();
  double initialAngle = getAngle();
  TinyVector<double,3> newPosMother, newPosDaughter;
  double separationDistance = 1.0*state_.cellLength0_/2.0;

  newPosMother(0) = initialPos(0) + cos(initialAngle*PI/180.0)*separationDistance;
  newPosMother(1) = initialPos(1) + sin(initialAngle*PI/180.0)*separationDistance;

  newPosDaughter(0) = initialPos(0) - cos(initialAngle*PI/180.0)*separationDistance;
  newPosDaughter(1) = initialPos(1) - sin(initialAngle*PI/180.0)*separationDistance;

  Cell::setPosition(newPosMother);
  daughterCell.Cell::setPosition(newPosDaughter);

  return daughterCell; // Return by value: makes a copy of daughterCell with the
                       // Cell copy constructor and return the copied Cell object.
}

//------------------------------------------------------------------------------

double Cell::getTimeDependentVolume(const double time) const
{
  return state_.getTimeDependentVolume(time);
}

//------------------------------------------------------------------------------









