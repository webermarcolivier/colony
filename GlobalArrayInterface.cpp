/***************************************************************************//**
 * Project: Colony
 *
 * \file    GlobalArrayInterface.cpp
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

#include "GlobalArrayInterface.h"

//------------------------------------------------------------------------------

GlobalArrayInterface::GlobalArrayInterface()
{}

//------------------------------------------------------------------------------

GlobalArrayInterface::~GlobalArrayInterface()
{}

//------------------------------------------------------------------------------

Array<double,1>& GlobalArrayInterface::getGlobalXConcCells(CellCollection* cellCollectionPtr)
{
  int nSpeciesMilieu = cellCollectionPtr->getMilieu().getNSpecies();
  globalXConcCells_.reference( globalXConc_.getGlobalArray() (blitz::Range( nSpeciesMilieu, blitz::toEnd )) );
  return globalXConcCells_;
}

//------------------------------------------------------------------------------

pair<int,int> GlobalArrayInterface::getLocalChannelAndCellIndices(int i) const
{
  pair<int,int> kj(globalPropensities_.getLocalIndicesFromGlobal(i));
  --(kj.first); // Change cell index range from 0,...,nCells to -1,...,nCells-1
                // -1 refers to the milieu, range 0,...,nCells-1 refers to the cells.
  return kj;
}

//------------------------------------------------------------------------------

void GlobalArrayInterface::getLocalArrayList(CellCollection* cellCollectionPtr)
{
  getLocalArrayListPropensities(cellCollectionPtr);
  getLocalArrayListX           (cellCollectionPtr);
  getLocalArrayListXConc       (cellCollectionPtr);
}

void GlobalArrayInterface::getLocalArrayListPropensities(CellCollection* cellCollectionPtr)
{
  // Erase the list of local arrays
  localArrayListPropensities_.resize(0);

  int i;

  LocalArrayTarget<double> target;

  target.ptr = & ((*cellCollectionPtr).getMilieu().getPropensities());

  target.size = (*cellCollectionPtr).getMilieu().getPropensities().size();

  localArrayListPropensities_.push_back(target);


  for (i=0; i<cellCollectionPtr->getNCells(); i++)
  {
    target.ptr = & ((*cellCollectionPtr)[i].getPropensities());

    target.size = (*cellCollectionPtr)[i].getPropensities().size();

    localArrayListPropensities_.push_back(target);
  }
}

void GlobalArrayInterface::getLocalArrayListX(CellCollection* cellCollectionPtr)
{
  // Erase the list of local arrays
  localArrayListX_.resize(0);

  int i;

  LocalArrayTarget<int> target;

  target.ptr = & ((*cellCollectionPtr).getMilieu().getX());

  target.size = (*cellCollectionPtr).getMilieu().getX().size();

  localArrayListX_.push_back(target);


  for (i=0; i<cellCollectionPtr->getNCells(); i++)
  {
    target.ptr = & ((*cellCollectionPtr)[i].getX());

    target.size = (*cellCollectionPtr)[i].getX().size();

    localArrayListX_.push_back(target);
  }
}

void GlobalArrayInterface::getLocalArrayListXConc(CellCollection* cellCollectionPtr)
{
  // Erase the list of local arrays
  localArrayListXConc_.resize(0);

  int i;

  LocalArrayTarget<double> target;

  target.ptr = & ((*cellCollectionPtr).getMilieu().getXConc());

  target.size = (*cellCollectionPtr).getMilieu().getXConc().size();

  localArrayListXConc_.push_back(target);


  for (i=0; i<cellCollectionPtr->getNCells(); i++)
  {
    target.ptr = & ((*cellCollectionPtr)[i].getXConc());

    target.size = (*cellCollectionPtr)[i].getXConc().size();

    localArrayListXConc_.push_back(target);
  }
}


//------------------------------------------------------------------------------

void GlobalArrayInterface::setLocalArraysReference
(
  CellCollection* cellCollectionPtr,
  list<Array<double,1> >& globalArrayPartsListPropensities,
  list<Array<int,1> >&    globalArrayPartsListX,
  list<Array<double,1> >& globalArrayPartsListXConc
)
{
  setLocalArraysReferencePropensities
  (
    cellCollectionPtr,
    globalArrayPartsListPropensities
  );
  setLocalArraysReferenceX
  (
    cellCollectionPtr,
    globalArrayPartsListX
  );
  setLocalArraysReferenceXConc
  (
    cellCollectionPtr,
    globalArrayPartsListXConc
  );
}

void GlobalArrayInterface::setLocalArraysReferencePropensities
(
  CellCollection* cellCollectionPtr,
  list<Array<double,1> >& globalArrayPartsListPropensities
)
{
  list<Array<double,1> >::iterator it =  globalArrayPartsListPropensities.begin();

  // Set reference of the propensities array of the milieu and advance to the
  // next element in the list.
  (*cellCollectionPtr).getMilieu().setReferencePropensities( (*it) );
  it++;

  // Set reference of the propensities array of the cells.
  int i;
  for (i = 0;
       it != globalArrayPartsListPropensities.end();
       it++, i++)
  {
    (*cellCollectionPtr)[i].setReferencePropensities( (*it) );
  }
}

void GlobalArrayInterface::setLocalArraysReferenceX
(
  CellCollection* cellCollectionPtr,
  list<Array<int,1> >& globalArrayPartsListX
)
{
  list<Array<int,1> >::iterator it =  globalArrayPartsListX.begin();

  // Set reference of the state array of the milieu and advance to the
  // next element in the list.
  (*cellCollectionPtr).getMilieu().setReferenceX( (*it) );
  it++;

  // Set reference of the state array of the cells.
  int i;
  for (i = 0;
       it != globalArrayPartsListX.end();
       it++, i++)
  {
    (*cellCollectionPtr)[i].setReferenceX( (*it) );
  }
}

void GlobalArrayInterface::setLocalArraysReferenceXConc
(
  CellCollection* cellCollectionPtr,
  list<Array<double,1> >& globalArrayPartsListXConc
)
{
    list<Array<double,1> >::iterator it =  globalArrayPartsListXConc.begin();

    // Set reference of the state array of the milieu and advance to the
    // next element in the list.
    (*cellCollectionPtr).getMilieu().setReferenceXConc( (*it) );
    it++;

    // Set reference of the state array of the cells.
    int i;
    for (i = 0;
         it != globalArrayPartsListXConc.end();
         it++, i++)
    {
      (*cellCollectionPtr)[i].setReferenceXConc( (*it) );
    }
}

//------------------------------------------------------------------------------

//void GlobalArrayInterface::buildArraysDirectMethod(CellCollection* cellCollectionPtr)
//{
//  getLocalArrayList(cellCollectionPtr);
//  globalPropensities_.build(localArrayList_);
//}

//------------------------------------------------------------------------------

void GlobalArrayInterface::buildArrays(CellCollection* cellCollectionPtr)
{
  /** - Get the list of local arrays. */
  getLocalArrayList(cellCollectionPtr);

  /** - Build the global array and get the list of global array ranges. */
  list<Array<double,1> > globalArrayRangesListPropensities =
  globalPropensities_.buildGlobalArrayRangesList(localArrayListPropensities_);

  list<Array<int,1> > globalArrayRangesListX =
  globalX_.buildGlobalArrayRangesList(localArrayListX_);

  list<Array<double,1> > globalArrayRangesListXConc =
  globalXConc_.buildGlobalArrayRangesList(localArrayListXConc_);

  /** - Call the setReferencePropensities in each cell of the cell collection
   *    with reference targets given in the previous list of global array ranges.
   */
  setLocalArraysReference(cellCollectionPtr,
                          globalArrayRangesListPropensities,
                          globalArrayRangesListX,
                          globalArrayRangesListXConc
                          );
}

//------------------------------------------------------------------------------








