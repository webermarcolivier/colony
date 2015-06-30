/***************************************************************************//**
 * Project: Colony
 *
 * \file    GlobalArrayInterface.h
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

#ifndef GLOBALARRAYINTERFACE_H
#define GLOBALARRAYINTERFACE_H

#include "debug.h"

// standard C++ header files
#include <iostream>
#include <list>
#include <utility>

// libraries header files
#include <blitz/array.h>

// user header files
#include "GlobalArray.h"

// namespaces
using blitz::Array;
using blitz::Range;
using std::cout;
using std::ostream;
using std::endl;
using std::list;
using std::pair;

/**
 * Use forwared declaration of class CellCollection. The header file
 * "CellCollection.h" is included
 * after the declaration of the GlobalArrayInterface class.
 */
class CellCollection;


/**
  * Interface between the class CellCollection and the class GlobalArray.
  *
  * This class allows the use of the GlobalArray class from the CellCollection
  * class. It provides with the methods needed to collect the local arrays
  * in each cell of the colony and to build a global array that regroup all
  * the local arrays. The local arrays are made to reference to a part
  * of the global array.
  */

class GlobalArrayInterface
{
public:

//---- LIFECYCLE

  GlobalArrayInterface();

  ~GlobalArrayInterface();

//---- OPERATORS

//---- ACCESS

  /**
   * Get the global array of propensities (by reference).
   */
  Array<double,1>& getGlobalPropensities();

  /**
   * Get the global array of x states (by reference).
   */
  Array<int,1>& getGlobalX();

  /**
   * Get the global array of xConc states (by reference).
   */
  Array<double,1>& getGlobalXConc();

  /**
   * Get the global array of xConc states only for the cells (by reference).
   */
  Array<double,1>& getGlobalXConcCells(CellCollection* cellCollectionPtr);

  /**
   * Get the cell index \c k and local channel index \c j
   * from the global channel index \c i.
   * @param i [in] global channel index.
   * @return indexPair [out] pair (cell index, local channel index).
   * The cell index \c k ranges from -1 to nCells-1. The index -1 refers to the
   * milieu and the index range 0,...,nCells-1 refers to the cells.
   * @see CellCollection#applyReaction()
   */
   pair<int,int> getLocalChannelAndCellIndices(int i) const;


//---- INQUIRY

//---- OPERATIONS

  /**
   * Construct the global arrays from the local ones.
   * Pass the list of local arrays to the GlobalArray member to build it.
   * @see GlobalArray
   */
//  void buildArraysDirectMethod(CellCollection* cellCollectionPtr);

  /**
   * Construct the global arrays from the local ones.
   * @see GlobalArray
   */
  void buildArrays(CellCollection* cellCollectionPtr);


private:

  // Please don't look at the private section, it is a mess...

//---- OPERATIONS

  /**
   * Get a list of all the local arrays in the cell collection that need to be
   * register in the global arrays. Writes the list in #localArrayList_.
   */
  void getLocalArrayList(CellCollection* cellCollectionPtr);
  void getLocalArrayListPropensities(CellCollection* cellCollectionPtr);
  void getLocalArrayListX(CellCollection* cellCollectionPtr);
  void getLocalArrayListXConc(CellCollection* cellCollectionPtr);

  /**
   * Set the reference of local arrays by calling the setReferencePropensities
   * method in the Cell objects.
   * @param [in,out] cellCollectionPtr Pointer to the cell collection.
   * @param [in] globalArrayPartsList List of the partial arrays of the global
   *        array (globalArray(Range(a,b)) used as reference targets for
   *        referencing the local arrays.
   */
  void setLocalArraysReference
  (
    CellCollection* cellCollectionPtr,
    list<Array<double,1> >& globalArrayPartsListPropensities,
    list<Array<int,1> >&    globalArrayPartsListX,
    list<Array<double,1> >& globalArrayPartsListXConc
  );

  // I warned you...

  void setLocalArraysReferencePropensities
  (
    CellCollection* cellCollectionPtr,
    list<Array<double,1> >& globalArrayPartsListPropensities
  );
  void setLocalArraysReferenceX
  (
    CellCollection* cellCollectionPtr,
    list<Array<int,1> >& globalArrayPartsListX
  );
  void setLocalArraysReferenceXConc
  (
    CellCollection* cellCollectionPtr,
    list<Array<double,1> >& globalArrayPartsListXConc
  );


//---- DATA

  /**
   * Global array of the propensities. It contains consecutively:
   * - propensities array of the milieu.
   * - propensities array of cells 0,...,nCells-1.
   *   - each cell array contains the propensities of internal
   *     reactions and cell-milieu reactions.
   */
  GlobalArray<double> globalPropensities_;

  /**
   * Global array of the states x. It contains consecutively:
   * - state array of the milieu.
   * - state array of cells 0,...,nCells-1.
   */
  GlobalArray<int> globalX_;
  GlobalArray<double> globalXConc_;

  /**
    * This array is just a reference to a part of globalXConc_.
    */
  Array<double,1> globalXConcCells_;

  /**
   * List of the local array targets (pointer+size).
   * This list is used to build the global array. The list form is just to
   * keep a general way of communication with the GlobalArray class.
   * @see LocalArrayTarget
   */
  list< LocalArrayTarget<double> > localArrayListPropensities_;
  list< LocalArrayTarget<int> >    localArrayListX_;
  list< LocalArrayTarget<double> > localArrayListXConc_;


};

//------------------------------------------------------------------------------

inline Array<double,1>& GlobalArrayInterface::getGlobalPropensities()
  {return globalPropensities_.getGlobalArray();}

inline Array<int,1>& GlobalArrayInterface::getGlobalX()
  {return globalX_.getGlobalArray();}

inline Array<double,1>& GlobalArrayInterface::getGlobalXConc()
  {return globalXConc_.getGlobalArray();}



#include "CellCollection.h"


#endif // GLOBALARRAYINTERFACE_H
