/***************************************************************************//**
 * Project: Colony
 *
 * \file    CellCollection.h
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

#ifndef CELLCOLLECTION_H
#define CELLCOLLECTION_H

#include "debug.h"

// standard C++ header files
#include <iostream>
#include <iomanip>
#include <utility>
#include <map>
#include <vector>
#include <string>
#include <limits> ///< numerical limits for each of the fundamental types


// libraries header files
#include <blitz/array.h>

// user header files

#include "Milieu.h"
#include "CellLineageGeneration.h"
#include "Cell.h"
#include "GlobalArrayInterface.h"

// namespaces
using blitz::Array;
using blitz::Range;
using std::cout;
using std::ostream;
using std::endl;
using std::setw;
using std::pair;
using std::multimap;
using std::vector;
using std::string;



/**
 * Use forward declaration of class Input. The header file "Input.h" is included
 * after the declaration of the CellCollection class.
 */
class Input;



class CellCollectionParam;



/**
 * A collection of cells.
 *
 * CellCollection contains a collection of cells and constitute the
 * framework for interacting with it and handling it. It defines the methods to
 * initialize a collection, access a particular cell by its index, delete and
 * duplicate one cell.
 *
 * It also defines an operation to compute the propensities of all reaction
 * channels in the collection and an operation to apply a reaction event.
*/
class CellCollection
{
public:


//---- LIFECYCLE

  CellCollection();

  ~CellCollection();

  /**
   * Initialize a homogeneous collection of cells. Each cell is initialized
   * with the same parameters. The milieu is also initialized.
   */
  void initializeCellCollection(CellCollectionParam& p);

  /**
   * Initialize a homogeneous collection of cells. Each cell is initialized
   * with the same parameters. The milieu is also initialized.
   * This version of the method takes as argument an Input object that contains
   * the parameters.
   */
  void initializeCellCollection(Input& input);


//---- OPERATORS

  /**
   * Operator that outputs the content of all cells.
   * Remark: it is declared as a friend operator, so that it has access to the
   * private part of the class.
   */
  friend ostream& operator<< (ostream& out, const CellCollection& cellCollection);

  /**
   * Random access to the cells.
   */
  const Cell& operator[] (int i) const;
        Cell& operator[] (int i);


//---- ACCESS

  /**
   * Non-constant access to the milieu.
   */
  Milieu& getMilieu();

  /**
   * Get the number of cells.
   */
  int getNCells() const;

  bool getConstantCellDensity() const;

  Array<double,1> getTimeDependentVolumes(const double time) const;

  double getTimeDependentVolumeMilieu(const double time) const;

  double updateTimeDependentVolumeMilieu(const double time);

  /**
   * Get the global array of propensities.
   *
   * Get the global array of propensities of the whole cell collection,
   * including the cell-milieu reaction channels.
   */
  Array<double,1>& getGlobalPropensities();

  Array<int,1>& getGlobalX();

  /**
   * Remark: XConc is not a global array, it is constructed from the concentration arrays of the cells.
   */
  Array<double,1> getXConc(const double time);

  /**
   * Get the global array of xConc. Note: this global array is only used when
   * we use the chemicalLangevin algorithm (then we describe the system with
   * species concentrations, not number of molecules).
   */
  Array<double,1>& getGlobalXConc();
  Array<double,1>& getGlobalXConcCells();

  double getSumPropensities() const;

  int getCellIndexFromGlobalChannelIndex(int mu);

  double getTimeNextDivision();

  Array<double,1> getCellCyclePhases(const double time) const;

  list<CellLineageGeneration> getCellLineage() const;

  Array<double,1> getCellsAngle() const;

  Array<TinyVector<double,3>,1> getCellsPosition() const;

  vector<string> getListSpeciesName() const;


//---- INQUIRY


//---- OPERATIONS


  void setCellsAngle(Array<double,1> cellsAngleArray);

  void setCellsPosition(Array<TinyVector<double,3>,1> cellsPositionArray);

  /**
   * Compute and get propensities of all reaction channels (including intercellular or
   * cell-milieu reactions).
   */
  Array<double,1>& computeAllCellsAndGetPropensities();

  /**
   * Compute and get propensities of all reaction channels (including intercellular or
   * cell-milieu reactions).
   * This method is optimized by using the #lastReactionCellIndex_.
   * Compute propensities ONLY in cell i: at each step,
   * we keep track of the cell index #lastReactionCellIndex_ where the reaction happens, and at next step,
   * we compute the propensities only in this cell and the milieu (because only
   * in these cells the internal state has changed).
   */
  Array<double,1>& computeAndGetPropensities();
  Array<double,1>& computeAndGetPropensitiesModified(const double time);

  void computeTimeDependentPropensities(const double time);
  void computePropensitiesAllInternalTimeDependentReactions(const double time);
  double getSumPropensitiesAllInternalTimeDependentReactions();

  /**
   * Initialize the sum of propensities.
   * We substract the sum of propensities from the milieu and set artificially
   * the lastReactionCellIndex to -1 (corresponding to the milieu),
   * in order to ensure the good initialization of the optimized sum of the
   * propensities.
   * This initialization has to be applied after dividing or deleting a cell.
   */
  void initializeSumPropensities();
  void initializeSumPropensities(const double time);

private:
  /**
   * Duplicate cell \c i.
   * The cell is duplicated and proteins are distributed between mother and
   * daughter cells with binomial probability.
   * - Substract the volume of the cell V0 to the volume of the milieu.
   *   (Remark: The dynamic variation of cells volume is not taken into account
   *   for keeping the milieu's volume constant.).
   * - Increase the size of the cells array with the resizeAndPreserve method.
   * - Call the duplicate method of cell \c i and copy the daughter cell to
   *   the new slot.
   * - nCells++
   * - Rebuild global arrays.
   * - Rebuild the list of the upcoming cell division events.
   */
  void duplicateCell(int i);

  /**
   * Compute the sum of the propensities of the cell-milieu reactions of
   * all cells.
   */
  double getSumPropensitiesAllCellsMilieuReactions();

  /**
   * Compute the propensities of the cell-milieu reactions of all cells.
   */
  void computePropensitiesCellMilieuAllCells();
  void computePropensitiesCellMilieuAllCells(const double time);


public:

  /**
   * Divide the cell(s) that are next in the division events list.
   * This method looks at the list of division events and divide the cell that
   * is next in the list. In the case that several cells have exactly the
   * same time of next division, we divide all these cells.
   */
  void applyNextDivisionEvent();

  /**
   * Delete cell \c i.
   * - Add the volume of the cell to the volume of the milieu.
   * - Make a copy of cells array.
   * - Resize the cells array to nCells-1.
   * - Copy back the cells (0,...,i-1,i+1,...,nCells) into the cells array.
   * - nCells--
   * - Update the cell generation list.
   * - Rebuild global arrays.
   * - Rebuild the list of the upcoming cell division events.
   */
  void deleteCell(int i, const double time);

  /**
   * Build the list of upcoming division events. This functions is called from
   * the #deleteCell() method. It iterates over the cells to get the value
   * of their timeNextDivision and their cell index and build the list.
   * It is necessary to rebuild the list when deleting one cell because almost
   * all the cell indices are changed after the deletion (by shifting).
   * Thus, it is simpler to rebuild the entire list than changing the value
   * of the cell index in each element.
   */
  void buildListDivisionEvents();

  /**
   * Apply a reaction corresponding to reaction channel \c mu (global index).
   * - Find cell number \c k and local channel number \c mu_k.
   * - If \c k = -1
   *   - Apply reaction \c mu_k in the milieu.
   * - else
   *   - Apply reaction \c mu_k in cell \c k.
   */
  void applyReaction(int mu);


private:

//---- DATA

  /**
   * Number of cells in the cell collection. Remark that this does not include
   * the milieu, which is considered as a special cell.
   */
  int nCells_;

  /**
   * The collection of cell objects is contained in a Blitz++ array. This will
   * ensure a good performance for random access to the cells (by index).
   * The drawback is that when deleting a cell, we have to perform complicated
   * copying operations. We suppose that duplicating/deleting events are rare,
   * so the lost in performance is negligable.
   */
  Array<Cell,1> cells_;

  /**
   * Milieu cell object that represents the milieu / environment where the cells
   * are living. It functions as a normal cell in the sense that it can define
   * some chemical species and reactions. Milieu can only have one instance.
   */
  Milieu milieu_;

  /**
   * The GlobalArrayInterface class defines the interface to the class that handle
   * global arrays. These are contiguous arrays that gather up
   * the internal arrays of all cells in order to be able to keep all cell data
   * in a unique place. This way, the global array can be passed easily to the
   * integrator interface and all the cell data is stored contiguously in memory.
   * Moreover, the local/internal arrays of the cells are referenced to a part of
   * the global array, i.e. they are sharing the same data. We have a
   * transparent manner of accessing the arrays, at the level of the cell
   * collection (global array) or at the level of the cell (local array).
   * Any changes made to the local array Cell#totalPropensities_ when calling
   * the Cell#computePropensities() method in a cell will also be reflected in the
   * global array, because they refer to the same data.
   *
   * Remark 1: When deleting or adding a cell, the global array need to be
   * completely rebuild.
   */
  GlobalArrayInterface globalArrayInterface_;

  /**
   * List of the upcoming cell division events. This list is used to know what
   * is the next cell divison event and when it will happen. It is implemented
   * using the multimap container from the STL. Each element in the multimap
   * has a key equal to the time of next division event with the cell index as
   * mapped value. The multimap always stays in a sorted order following the
   * < comparison for the key values (times of division).
   *
   * Remark: The cell index refers to the cells array in the state when the list
   * was built. When applying a division event, the array has one more element
   * and its size is increased to nCells+1. However, the new daughter cell is
   * inserted at the very end of the array, thus indices from 0,...,nCells-1
   * can still be used to refer to the cells and are consistent with the indices
   * before the division. This is especially important when two or more cells
   * have the same division event time and are divided at the same time.
   */
  multimap<double,int> listDivisionEvents_;

  /**
   * List that contains information about the cell lineage.
   * @see CellLineageGeneration.
   */
  list<CellLineageGeneration> cellLineage_;

  /**
   * Temporary variable: index of the cell in which the reaction of previous integration step
   * happened. This index is used to optimize the integration algorithm:
   * @see #computeAndGetPropensities()
   * @see #computeSumPropensities()
   *
   * Remark: if the reaction happened in the milieu, it is set to -1.
   */
  int lastReactionCellIndex_;

  /**
   * Temporary variable: indicates if the reaction that happened
   * at the previous integration step was an internal (cell) reaction or an
   * external (cell-milieu) reaction. This index is used to optimize the integration algorithm:
   * @see #computeAndGetPropensities()
   * @see #computeSumPropensities()
   *
   * Remark: if the reaction happened in the milieu, then it cannot be external
   * (by definition, the milieu object does not have any external reactions).
   */
  bool isLastReactionInternal_;

//  /**
//   * Temporary variable: true if the last reaction was an intracellular and
//   * false if it was a cell-milieu reaction.
//   */
//  bool lastReactionIsIntracellular_;

  /**
   * Temporary variable: sum of the propensities. Used in the optimization
   * of the computation of the sum of all the propensities. The idea of the
   * optimization stands as follows:
   * - The next reaction to occur introduces
   *   a change in the state of the system, but only in the number of molecules
   *   of the species involved in the reaction.
   *   We assume that there is only a reduced set of reaction channels
   *   whose propensities depend on the
   *   species number that has changed during the reaction. Then, we only need to
   *   update the propensities of the reaction channels that are affected by
   *   the last reaction.
   *
   * - The propensities and the sum are initialized. Then, in the integration loop:
   * - Before applying the reaction:
   *   - Propensities of reaction channels that are affected by the next reaction
   *     are substracted from the sum.
   *   - Apply reaction.
   * - Compute propensities:
   *   - Compute propensities of the channels affected by the last reaction.
   *   - Add the sum of the propensities of these channels to the sum.
   * - The propensities and their sum have been updated.
   *
   * In this code, we only use the connections of the reactions between
   * cell-cell and cell-milieu. We do not analyse the internal chemical system in the
   * cell and the connections between reactions.
   */
  double sumPropensities_;

  /**
   * If true the cell density is kept constant by automatically killing the
   * cells (at random) when cells duplicate.
   */
  bool constantCellDensity_;

};

//------------------------------------------------------------------------------

// INLINE ACCESS METHODS

inline int CellCollection::getNCells() const
  {return nCells_;}

inline bool CellCollection::getConstantCellDensity() const
  {return constantCellDensity_;}

inline Array<double,1>& CellCollection::getGlobalPropensities()
  {return globalArrayInterface_.getGlobalPropensities();}

inline Array<int,1>& CellCollection::getGlobalX()
  {return globalArrayInterface_.getGlobalX();}

inline Array<double,1>& CellCollection::getGlobalXConc()
  {return globalArrayInterface_.getGlobalXConc();}

inline Array<double,1>& CellCollection::getGlobalXConcCells()
  {
    return globalArrayInterface_.getGlobalXConcCells(this);
  }

inline double CellCollection::getSumPropensities() const
  {return sumPropensities_;}

inline Milieu& CellCollection::getMilieu()
  {return milieu_;}

inline double CellCollection::getTimeNextDivision()
{
  #ifdef TIME_DEPENDENT_PROPENSITIES
    multimap<double,int>::iterator it (listDivisionEvents_.begin());
    pair<double,int> nextDivisionEvent (*it);
    return nextDivisionEvent.first;
  #else
    return numeric_limits<double>::infinity();
  #endif
}

inline list<CellLineageGeneration> CellCollection::getCellLineage() const
  {return cellLineage_;}


#include "Param/CellCollectionParam.h"
#include "Input.h"


#endif // CELLCOLLECTION_H
