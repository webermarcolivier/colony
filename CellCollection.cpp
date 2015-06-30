/***************************************************************************//**
 * Project: Colony
 *
 * \file    CellCollection.cpp
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

#include "CellCollection.h"

//------------------------------------------------------------------------------

CellCollection::CellCollection()
  : nCells_(0),
    lastReactionCellIndex_(-1),
    isLastReactionInternal_(true)
{}

//------------------------------------------------------------------------------

CellCollection::~CellCollection()
{
  cells_.free();
}

//------------------------------------------------------------------------------

ostream& operator<< (ostream& out, const CellCollection& cellCollection)
{
  int i;
  int nCells = cellCollection.getNCells();
  out << "################\n";
  out << "Cell collection:\n";
  out << "Number of cells in cell collection = " << nCells << '\n';
  out << "Total number of Cell objects = " << Cell::getNCellObjects() << "\n\n";
  out << "Milieu:\n" << cellCollection.milieu_;

  for (i=0; i<nCells; i++)
  {
    out << "Cell#" << setw(4) << i << ":\n" << cellCollection[i];
  }
  out << "################\n";
  return out;
}

//------------------------------------------------------------------------------

const Cell& CellCollection::operator[] (int i) const
{
  return cells_(i);
}

      Cell& CellCollection::operator[] (int i)
{
  return cells_(i);
}

//------------------------------------------------------------------------------

void CellCollection::initializeCellCollection(CellCollectionParam& p)
{

  // Calculate volume of the milieu from the total volume and the volume of
  // the cells.
  double volume0Cell = p.cellInitParam_.cellBaseInitParam_.stateInitParam_.cellVolume0_;

  // UPDATE: Now the volume0 for the milieu is the total volume and we substract always
  // the sum of the volume of the cells to the total volume to get the volume of the milieu
  //double milieuVolume0 = p.totalVolume_ - p.nCells_*volume0Cell;
  double milieuVolume0 = p.totalVolume_;

  // Initialize the milieu.
  p.milieuInitParam_.stateInitParam_.cellVolume0_ = milieuVolume0;
  milieu_.initialize(p.milieuInitParam_);
  updateTimeDependentVolumeMilieu(0.0);

  // Set the cells to point to the milieu.
  p.cellInitParam_.cellMilieuChemicalSystemInitParam_.milieuPtr_ = &milieu_;

  // Initialize the array of cells.
  int nCells = p.nCells_;
  cells_.resize(nCells);
  nCells_ = nCells;


  if (!p.x0CellHeterogenous_)
  {
    // Initialize cell collection with the same initial conditions for all cells.

    int i;
    for (i=0; i<nCells; i++)
    {
      #ifdef RANDOM_INITIAL_PHASES
        // Setting a random distribution of cell cycle initial phases.
        double phase = RandomNumberGenerator::getUniform();
        p.cellInitParam_.cellCyclephase_ = phase;
      #else
        // Setting the cell cycle initial phases for all cells.
        // We set the initial phase such that the initial volume of
        // the cells is equal to V1, the cell volume averaged over
        // a cell cycle, V1 = V0 / ln(2). The exact value of the
        // inial phase is ph1 = - ln( ln(2) ) / ln(2).
        p.cellInitParam_.cellCyclephase_ = 0.528766372944898;
      #endif

      cells_(i).initialize(p.cellInitParam_);
      cells_(i).setCellIndex(i);
    }

  } else {
    // Initialize cell collection with different initial conditions for every cell.

    int i;
    for (i=0; i<nCells; ++i)
    {
      // Copy the initial conditions (cell cycle phase and number of molecules of chemical species)
      // in the CellInitParam class
      p.cellInitParam_.cellCyclephase_ = p.initialPhaseTable_(i);
      p.cellInitParam_.cellBaseInitParam_.stateInitParam_.x_ = p.x0CellHeterogenousTable_(i, blitz::Range::all());
      // Initialize the cell #i with the copied values
      cells_(i).initialize(p.cellInitParam_);
      cells_(i).setCellIndex(i);
    }
  }


  // Initialize the global array of propensities.
  globalArrayInterface_.buildArrays(this);

  // Initialize the list of upcoming division events.
  buildListDivisionEvents();

  #ifndef USE_CHEMICAL_LANGEVIN
    // Initialize the temporary sum of propensities.
    #ifdef TIME_DEPENDENT_PROPENSITIES
      initializeSumPropensities(0.0);
    #else
      initializeSumPropensities();
    #endif // TIME_DEPENDENT_PROPENSITIES
  #endif // USE_CHEMICAL_LANGEVIN


  // Initialize the cell lineage list.
  cellLineage_.clear();
  CellLineageGeneration cellLineageGeneration;
  Array<int,1> motherCellsIndices(nCells);
  int i;
  for (i=0; i<nCells; i++)
  {
    motherCellsIndices(i) = i;
  }
  cellLineageGeneration.initialize(0.0, nCells, 0, motherCellsIndices);
  CellLineageGeneration cellLineageGenerationCopy( cellLineageGeneration );
  cellLineage_.push_back(cellLineageGenerationCopy);

  // Initialize the boolean parameter
  constantCellDensity_ = p.constantCellDensity_;

  // Initialize cells positions and angle.
  TinyVector<double,3> position;
  double cellLength0 = p.cellInitParam_.cellBaseInitParam_.stateInitParam_.cellLength0_;
  double xStep = 1.5*cellLength0;
  int k, l;
  int xn = int(sqrt(float(nCells_)));
  double x, y;
  for(i=0;i<nCells;++i)
  {
    // Set the same angle for all cells.
    (*this)[i].setAngle( 35.0 );
    // Regular distribution
    k = i / xn;
    l = i%xn;
    x = (l+1-xn/2)*xStep;
    y = (k+1-xn/2)*xStep;
    position(0) = x;
    position(1) = y;
    position(2) = 0;
    (*this)[i].setPosition(position);
    // Update the graphics cell size and position
    (*this)[i].updateGraphicsCell(0.0);
  }

}

//------------------------------------------------------------------------------

void CellCollection::initializeCellCollection(Input& input)
{
  initializeCellCollection(input.p);
}

//------------------------------------------------------------------------------

Array<double,1>& CellCollection::computeAndGetPropensities()
{
  #ifdef OPTIMIZE_COMPUTE_PROPENSITIES

    if ( lastReactionCellIndex_ != -1 )
    {
      // Last reaction happened in cell k.

      if (isLastReactionInternal_)
      {
        // Compute the propensities of internal and external reactions of cell k.
        cells_(lastReactionCellIndex_).computePropensities();

        sumPropensities_ += sum( (*this)[lastReactionCellIndex_].getPropensities() );
      }
      else
      {
        // Compute the propensities of internal reactions of cell k,
        // of the milieu reactions and the external reactions of all cells.
        cells_(lastReactionCellIndex_).computePropensitiesCell();

        computePropensitiesCellMilieuAllCells();

        milieu_.computePropensities();

        sumPropensities_ += sum( (*this)[lastReactionCellIndex_].getPropensitiesCell() ) +
                            sum( milieu_.getPropensities() ) +
                            getSumPropensitiesAllCellsMilieuReactions();
      }

    } else
    {
      // Last reaction happened in the milieu.

      milieu_.computePropensities();

      computePropensitiesCellMilieuAllCells();

      sumPropensities_ += sum( milieu_.getPropensities() ) +
                          getSumPropensitiesAllCellsMilieuReactions();
    }

  #else

    int i;

    milieu_.computePropensities();

    for (i=0; i<nCells_; i++)
    {
      cells_(i).computePropensities();
    }

    sumPropensities_ = sum( getGlobalPropensities() );

  #endif // OPTIMIZE_COMPUTE_PROPENSITIES


  return getGlobalPropensities();
}
//------------------------------------------------------------------------------

Array<double,1>& CellCollection::computeAndGetPropensitiesModified
(
  const double time
)
{
  #ifdef OPTIMIZE_COMPUTE_PROPENSITIES

// TODO (mweber#1#): Try to figure out how to build an optimization with time-dependent reaction channels everywhere (in cells, Milieu, and cell-milieu reactions).

    if ( lastReactionCellIndex_ != -1 )
    {
      // Last reaction happened in cell k.

      if (isLastReactionInternal_)
      {
        // Compute the propensities of internal and external reactions of cell k.
        cells_(lastReactionCellIndex_).computePropensities(time);

        sumPropensities_ += sum( (*this)[lastReactionCellIndex_].getPropensities() );
      }
      else
      {
        // Compute the propensities of internal reactions of cell k and
        // the external reactions of all cells.
        cells_(lastReactionCellIndex_).computePropensitiesCell(time);

        computePropensitiesCellMilieuAllCells(time);

        milieu_.computePropensities(time);

        sumPropensities_ += sum( (*this)[lastReactionCellIndex_].getPropensitiesCell() ) +
                            sum( milieu_.getPropensities() ) +
                            getSumPropensitiesAllCellsMilieuReactions();
      }

    } else
    {
      // Last reaction happened in the milieu.

      milieu_.computePropensities(time);

      computePropensitiesCellMilieuAllCells(time);

      sumPropensities_ += sum( milieu_.getPropensities() ) +
                          getSumPropensitiesAllCellsMilieuReactions();
    }

  #else

    int i;

    milieu_.computePropensities(time);

    for (i=0; i<nCells_; i++)
    {
      cells_(i).computePropensities(time);
    }

    sumPropensities_ = sum( getGlobalPropensities() );

  #endif // OPTIMIZE_COMPUTE_PROPENSITIES

  return getGlobalPropensities();
}

//------------------------------------------------------------------------------

void CellCollection::computePropensitiesAllInternalTimeDependentReactions(const double time)
{
  // TO IMPLEMENT (SEE CELLBASE)
  exit(1);
  int i;
  for (i=0; i<nCells_; i++)
  {
    (*this)[i].computePropensitiesTimeDependentReactions(time);
  }
}

//------------------------------------------------------------------------------


void CellCollection::computeTimeDependentPropensities(const double time)
{
  // TO IMPLEMENT (SEE CELLBASE)
  exit(1);
  #ifdef OPTIMIZE_COMPUTE_PROPENSITIES

    sumPropensities_ -= getSumPropensitiesAllCellsMilieuReactions();
    sumPropensities_ -= getSumPropensitiesAllInternalTimeDependentReactions();

    computePropensitiesCellMilieuAllCells(time);
    computePropensitiesAllInternalTimeDependentReactions(time);

    sumPropensities_ += getSumPropensitiesAllCellsMilieuReactions();
    sumPropensities_ += getSumPropensitiesAllInternalTimeDependentReactions();

  #else

    computePropensitiesCellMilieuAllCells(time);
    computePropensitiesAllInternalTimeDependentReactions(time);

    sumPropensities_ = sum( getGlobalPropensities() );

  #endif // OPTIMIZE_COMPUTE_PROPENSITIES
}

//------------------------------------------------------------------------------

void CellCollection::duplicateCell(int i)
{
  //double volume0 = (*this)[i].getVolume0();
  //milieu_.changeVolume0By(-volume0); //UPDATE: now volume0 for milieu is constant

  cells_.resizeAndPreserve(nCells_+1); // IMPORTANT: preserve the cells 0,..,nCells-1 !!!

  cells_(nCells_) = (*this)[i].duplicate();
  cells_(nCells_).setCellIndex(nCells_);

  nCells_++;
}

//------------------------------------------------------------------------------

void CellCollection::applyNextDivisionEvent()
{
  /*
   * Remark: the cell indices are still consistent after a division event
   * because the daughter cells are inserted at the end of the cell array.
   */

  int nCellsBefore = getNCells();

  // Point the iterator to the first division event in the list.
  multimap<double,int>::iterator it = listDivisionEvents_.begin();

  // Get the timeNextDivison of the first element.
  double timeNextDivision = (*it).first;

  // Declare a pair of iterators.
  pair<multimap<double,int>::iterator,multimap<double,int>::iterator>
    itRange;

  /* Returns the bounds of the range that includes all the elements in the
     container with a key that compares equal to x.
     If x does not match any key in the container, the range returned has a
     length of zero, with both iterators pointing to the element with nearest
     key greater than x, if any, or to multimap::end if x is greater than all
     the elements in the container. */
  itRange = listDivisionEvents_.equal_range( timeNextDivision );

  #ifdef BUILD_DEBUG
    // In principle, the following should never happen, since there is always
    // one element in the list and we perform a search in the list for the key
    // of the first element.
    if (itRange.first == itRange.second)
    {
      cout << "ERROR: class CellCollection, function applyNextDivisionEvent(),"
              " the iterator range of 'listDivisionEvents_' could not find the "
              "range of cells that have the same timeNextDivision." << endl;
      exit(1);
    }
  #endif // BUILD_DEBUG

  // Divide all cells that have the same time of next division as the first one
  // in the list.
  for (it=itRange.first; it!=itRange.second; ++it)
  {
    duplicateCell( (*it).second );
  }

  #ifdef PRINT_CELL_DIVISION_EVENT
    cout << "cell division event: " << "nCells = " << setw(4) << nCells_ <<
    " time = " << timeNextDivision << endl;
  #endif // PRINT_CELL_DIVISION_EVENT


  // Insert the cell lineage generation in the list.
  int nCellsAfter = getNCells();
  CellLineageGeneration cellLineageGeneration;
  Array<int,1> motherCellsIndices(nCellsAfter);
  // The new cells are added at the end of the array, thus the first part of the
  // indices are conserved.
  int i;
  for (i=0; i<nCellsBefore; i++)
  {
    motherCellsIndices(i) = i;
  }
  // The new indices are ordered as in the list of division events.
  for (it=itRange.first, i=nCellsBefore; it!=itRange.second; ++it, ++i)
  {
    motherCellsIndices(i) = (*it).second;
  }
  cellLineageGeneration.initialize
  (
    timeNextDivision,
    nCellsAfter,
    nCellsAfter-nCellsBefore,
    motherCellsIndices
  );
  CellLineageGeneration cellLineageGenerationCopy( cellLineageGeneration );
  cellLineage_.push_back(cellLineageGenerationCopy);


  // Rebuild the global array of propensities.
  globalArrayInterface_.buildArrays(this);

  // Rebuild the list of upcoming division events.
  buildListDivisionEvents();

  // Initialize the temporary sum of propensities.
  #ifdef TIME_DEPENDENT_PROPENSITIES
    initializeSumPropensities(timeNextDivision);
  #else
    initializeSumPropensities();
  #endif // TIME_DEPENDENT_PROPENSITIES

  // If we keep the cell density constant, kill the same number of cells that
  // has duplicated. Choose the cells at random.
  if (constantCellDensity_)
  {
    int nNewCells = nCellsAfter-nCellsBefore;
    int j;
    for (j=0;j<nNewCells;j++)
    {
      #ifdef RULE_DELETE_FARTHEST_CELL
        // Kill the most outer cell. Remark: This is **WRONG** statistically because the killed
        // cells are related to each other by cell division. We test it just for fun and
        // for making nicer movies of the simulations.
        int iCellOuter = 0;
        double distanceOuter = 0.0;
        for (int ik=0; ik<nCells_; ik++)
        {
          // Measuring the distance from the center of the colony (0,0).
          double d2 = (*this)[ik].getPosition()(0)*(*this)[ik].getPosition()(0)
                       + (*this)[ik].getPosition()(1)*(*this)[ik].getPosition()(1);
          if ( d2 > distanceOuter )
          {
            iCellOuter = ik;
            distanceOuter = d2;
          }
        }
        int iKill = iCellOuter;
      #else
        int iKill = floor( RandomNumberGenerator::getUniform() * getNCells() );
      #endif
      deleteCell(iKill, timeNextDivision);
    }
  }

  updateTimeDependentVolumeMilieu(timeNextDivision);
}

//------------------------------------------------------------------------------

void CellCollection::deleteCell(int iCell, const double time)
{
  //cout << "CellCollection::deleteCell " << iCell << endl;
  int nCellsBefore = getNCells();

  double volume0 = (*this)[iCell].getVolume0();

  milieu_.changeVolume0By(+volume0);

  Array<Cell,1> cellsCopy( cells_.copy() );

  cells_.resize(nCells_ - 1);

  int j;
  for (j=0; j<iCell; j++)
  {
    cells_(j) = cellsCopy(j);
  }
  for (j=iCell+1; j<nCells_; j++)
  {
    cells_(j-1) = cellsCopy(j);
  }

  for (j=0; j<nCells_-1; j++)
  {
    cells_(j).setCellIndex(j);
  }

  nCells_--;


  // Insert the cell lineage generation in the list.
  int nCellsAfter = getNCells();
  CellLineageGeneration cellLineageGeneration;
  Array<int,1> motherCellsIndices(nCellsAfter);
  // When the cell j is deleted, the array of cells is shifted for cells
  // j+1,...,nCells with a shift -1.
  for (j=0; j<iCell; j++)
  {
    motherCellsIndices(j) = j;
  }
  for (j=iCell; j<nCells_; j++)
  {
    motherCellsIndices(j) = j+1;
  }
  cellLineageGeneration.initialize
  (
    time,
    nCellsAfter,
    nCellsAfter-nCellsBefore,
    motherCellsIndices
  );
  CellLineageGeneration cellLineageGenerationCopy( cellLineageGeneration );
  cellLineage_.push_back(cellLineageGenerationCopy);


  // Rebuild the global array of propensities.
  globalArrayInterface_.buildArrays(this);

  // Rebuild the list of the upcoming cell division events.
  buildListDivisionEvents();

  // Initialize the temporary sum of propensities.
  #ifdef TIME_DEPENDENT_PROPENSITIES
    initializeSumPropensities(time);
  #else
    initializeSumPropensities();
  #endif


  updateTimeDependentVolumeMilieu(time);
}

//------------------------------------------------------------------------------

void CellCollection::buildListDivisionEvents()
{
  listDivisionEvents_.clear();

  multimap<double,int>::iterator it;

  int i;
  double timeNextDivision;

  for (i=0; i<nCells_; i++)
  {
    timeNextDivision = (*this)[i].getTimeNextDivision();

    it = listDivisionEvents_.insert( pair<double,int>(timeNextDivision,i) );
  }
}

//------------------------------------------------------------------------------

void CellCollection::applyReaction(int mu)
{
  // Find cell number k and local reaction channel mu_k.
  pair<int,int> indexPair
    ( globalArrayInterface_.getLocalChannelAndCellIndices(mu) );

  int k = indexPair.first;
  int mu_k = indexPair.second;

  #ifdef OPTIMIZE_COMPUTE_PROPENSITIES
    // Keep trace of the cell index where the reaction happened.
    lastReactionCellIndex_ = k;
  #endif //OPTIMIZE_COMPUTE_PROPENSITIES

  if (k != -1)
  {
    #ifdef OPTIMIZE_COMPUTE_PROPENSITIES
      // Keep trace if the reaction was internal or external.
      isLastReactionInternal_ = (*this)[k].isReactionInternal(mu_k);

      // Substract the sum of the propensities of the reactions channels affected
      // by the change of state.
      if (isLastReactionInternal_)
      {
        sumPropensities_ -= sum( (*this)[k].getPropensities() );
      } else
      {
        sumPropensities_ -= sum( (*this)[k].getPropensitiesCell() ) +
                            sum( milieu_.getPropensities() ) +
                            getSumPropensitiesAllCellsMilieuReactions();
      }
    #endif // OPTIMIZE_COMPUTE_PROPENSITIES

    // Apply reaction mu_k in cell k.
    (*this)[k].applyReaction(mu_k);

  } else
  {

    #ifdef OPTIMIZE_COMPUTE_PROPENSITIES
      // Keep trace if the reaction was internal or external.
      isLastReactionInternal_ = true;

      // Substract the sum of the propensities of the reactions channels affected
      // by the change of state.
      sumPropensities_ -= sum( milieu_.getPropensities() ) +
                          getSumPropensitiesAllCellsMilieuReactions();
    #endif // OPTIMIZE_COMPUTE_PROPENSITIES

    // Apply reaction mu_k in the milieu.
    milieu_.applyReaction(mu_k);
  }

}

//------------------------------------------------------------------------------

int CellCollection::getCellIndexFromGlobalChannelIndex(int mu)
{
  // Find cell number k.
  pair<int,int> indexPair
    ( globalArrayInterface_.getLocalChannelAndCellIndices(mu) );

  int k = indexPair.first;
  return k;
}

//------------------------------------------------------------------------------

void CellCollection::initializeSumPropensities()
{
  // Remark: we need to have already rebuilt the global array of propensities
  // before calling this method.

  // Compute all the propensities.
  int i;
  milieu_.computePropensities();
  for (i=0; i<nCells_; i++)
  {
    cells_(i).computePropensities();
  }

  // Sum the propensities and substract the propensities corresponding to the
  // artificial milieu reaction.
  sumPropensities_ = sum( getGlobalPropensities() )
                     - sum( milieu_.getPropensities() )
                     - getSumPropensitiesAllCellsMilieuReactions();
  lastReactionCellIndex_ = -1;
  isLastReactionInternal_ = true;
}

//------------------------------------------------------------------------------

void CellCollection::initializeSumPropensities(const double time)
{
  // Remark: we need to have already rebuilt the global array of propensities
  // before calling this method.

  // Compute all the propensities.
  int i;
  milieu_.computePropensities(time);
  for (i=0; i<nCells_; i++)
  {
    cells_(i).computePropensities(time);
  }

  // Sum the propensities and substract the propensities corresponding to the
  // artificial milieu reaction.
  sumPropensities_ = sum( getGlobalPropensities() )
                     - sum( milieu_.getPropensities() )
                     - getSumPropensitiesAllCellsMilieuReactions();
  lastReactionCellIndex_ = -1;
  isLastReactionInternal_ = true;
}

//------------------------------------------------------------------------------

double CellCollection::getSumPropensitiesAllCellsMilieuReactions()
{
  int i;
  double sumP = 0.0;
  for (i=0; i<nCells_; i++)
  {
    sumP += sum( (*this)[i].getPropensitiesCellMilieu() );
  }
  return sumP;
}

//------------------------------------------------------------------------------

void CellCollection::computePropensitiesCellMilieuAllCells()
{
  int i;
  for (i=0; i<nCells_; i++)
  {
    (*this)[i].computePropensitiesCellMilieu();
  }
}

//------------------------------------------------------------------------------

void CellCollection::computePropensitiesCellMilieuAllCells(const double time)
{
  int i;
  for (i=0; i<nCells_; i++)
  {
    (*this)[i].computePropensitiesCellMilieu(time);
  }
}

//------------------------------------------------------------------------------

Array<double,1> CellCollection::getCellCyclePhases(const double time) const
{
  Array<double,1> cellCyclePhases(nCells_);
  int i;
  for (i=0; i<nCells_; i++)
  {
    cellCyclePhases(i) = (*this)[i].getCellCyclePhase(time);
  }
  return cellCyclePhases;
}

//------------------------------------------------------------------------------

Array<double,1> CellCollection::getTimeDependentVolumes(const double time) const
{
  Array<double,1> volumeArray(nCells_);
  int i;
  for (i=0; i<nCells_; ++i)
  {
    volumeArray(i) = (*this)[i].getTimeDependentVolume(time);
  }
  return volumeArray;
}

//------------------------------------------------------------------------------

double CellCollection::getTimeDependentVolumeMilieu(const double time) const
{
  //return milieu_.getVolume0() - sum(getTimeDependentVolumes(time));
    return milieu_.getTimeDependentVolume(time);
}

//------------------------------------------------------------------------------

double CellCollection::updateTimeDependentVolumeMilieu(const double time)
{
  milieu_.milieuVolume_ = milieu_.getVolume0() - sum(getTimeDependentVolumes(time));
}

//------------------------------------------------------------------------------

Array<double,1> CellCollection::getXConc(const double time)
{
  Array<double,1> xConc;
  int nSpecies;
  nSpecies = milieu_.getNSpecies();
  xConc.resize(nSpecies);
  xConc(Range(0,nSpecies-1)) = milieu_.getXconc(time);

  int j0 = nSpecies;
  int i;
  for (i=0; i<nCells_; ++i)
  {
    nSpecies = (*this)[i].getNSpecies();
    xConc.resizeAndPreserve(xConc.size()+nSpecies);
    xConc(Range(j0,j0 + nSpecies - 1)) = (*this)[i].getXconc(time);
    j0 += nSpecies;
  }

  return xConc;
}

//------------------------------------------------------------------------------

Array<double,1> CellCollection::getCellsAngle() const
{
  Array<double,1> cellsAngleArray;
  cellsAngleArray.resize(nCells_);
  int i;
  for(i=0;i<nCells_;++i)
  {
    cellsAngleArray(i) = (*this)[i].getAngle();
  }
  return cellsAngleArray;
}

//------------------------------------------------------------------------------

Array<TinyVector<double,3>,1> CellCollection::getCellsPosition() const
{
  Array<TinyVector<double,3>,1> cellsPositionArray;
  cellsPositionArray.resize(nCells_);
  int i;
  for(i=0;i<nCells_;++i)
  {
    cellsPositionArray(i) = (*this)[i].CellBase::getPosition();
  }
  return cellsPositionArray;
}

//------------------------------------------------------------------------------

void CellCollection::setCellsAngle(Array<double,1> cellsAngleArray)
{
  if (cellsAngleArray.extent(firstDim) != nCells_)
  {
    cout << "ERROR: class CellCollection: method setCellsAngle, size of cellsAngleArray is not the same as nCells_." << endl;
    exit(1);
  }

  int i;
  for(i=0;i<nCells_;++i)
  {
    (*this)[i].setAngle(cellsAngleArray(i));
  }
}

//------------------------------------------------------------------------------

void CellCollection::setCellsPosition(Array<TinyVector<double,3>,1> cellsPositionArray)
{
  if (cellsPositionArray.extent(firstDim) != nCells_)
  {
    cout << "ERROR: class CellCollection: method setCellsPosition, size of cellsPositionArray is not the same as nCells_." << endl;
    exit(1);
  }

  int i;
  for(i=0;i<nCells_;++i)
  {
    (*this)[i].setPosition(cellsPositionArray(i));
  }
}

//------------------------------------------------------------------------------

vector<string> CellCollection::getListSpeciesName() const
{
  vector<string> listSpeciesName = (*this)[0].getListSpeciesName();
  return listSpeciesName;
}

//------------------------------------------------------------------------------
